// SPDX-FileCopyrightText: 2026 ficool2
// SPDX-License-Identifier: MIT

#include "fluid.h"

#include "body.h"
#include "core.h"
#include "physics_world.h"
#include "recording.h"
#include "sensor.h"
#include "shape.h"
#include "solver_set.h"

#include "arena_allocator.h"
#include "box3d/box3d.h"
#include "box3d/collision.h"

#include <float.h>
#include <stdlib.h>

//
// Based on "Exact Buoyancy for Polyhedra" by Erin Catto in Game Programming Gems 6 ;)
//

#define B3_FLUID_MAX_PLANES 32
#define B3_FLUID_MAX_POINTS 1024
#define B3_FLUID_MAX_FACES 256

typedef struct b3ClipPolyhedron
{
	b3Vec3 points[B3_FLUID_MAX_POINTS];
	int faceStarts[B3_FLUID_MAX_FACES];
	int faceCounts[B3_FLUID_MAX_FACES];
	int pointCount;
	int faceCount;
} b3ClipPolyhedron;

static void b3BeginFace( b3ClipPolyhedron* poly )
{
	if ( poly->faceCount < B3_FLUID_MAX_FACES )
	{
		poly->faceStarts[poly->faceCount] = poly->pointCount;
		poly->faceCounts[poly->faceCount] = 0;
	}
}

static bool b3PushPoint( b3ClipPolyhedron* poly, b3Vec3 p )
{
	if ( poly->pointCount == B3_FLUID_MAX_POINTS || poly->faceCount == B3_FLUID_MAX_FACES )
	{
		return false;
	}

	poly->points[poly->pointCount++] = p;
	poly->faceCounts[poly->faceCount] += 1;
	return true;
}

static void b3EndFace( b3ClipPolyhedron* poly )
{
	if ( poly->faceCount == B3_FLUID_MAX_FACES )
	{
		return;
	}

	if ( poly->faceCounts[poly->faceCount] < 3 )
	{
		// degenerate
		poly->pointCount = poly->faceStarts[poly->faceCount];
		return;
	}

	poly->faceCount += 1;
}

static bool b3MakeHullPolyhedron( b3ClipPolyhedron* poly, const b3HullData* hull, b3Transform transform )
{
	const b3HullHalfEdge* edges = b3GetHullEdges( hull );
	const b3HullFace* faces = b3GetHullFaces( hull );
	const b3Vec3* points = b3GetHullPoints( hull );

	poly->pointCount = 0;
	poly->faceCount = 0;

	for ( int i = 0; i < hull->faceCount; ++i )
	{
		b3BeginFace( poly );

		int first = faces[i].edge;
		int edge = first;
		int guard = 0;
		do
		{
			if ( b3PushPoint( poly, b3TransformPoint( transform, points[edges[edge].origin] ) ) == false )
			{
				return false;
			}

			edge = edges[edge].next;
			guard += 1;
		}
		while ( edge != first && guard < hull->edgeCount );

		b3EndFace( poly );
	}

	return poly->faceCount >= 4;
}

typedef struct b3CapSort
{
	float angle;
	int index;
} b3CapSort;

static int b3CompareCapAngles( const void* a, const void* b )
{
	const b3CapSort* ca = (const b3CapSort*)a;
	const b3CapSort* cb = (const b3CapSort*)b;
	return ( ca->angle > cb->angle ) - ( ca->angle < cb->angle );
}

typedef struct b3ClipBuffers
{
	b3ClipPolyhedron polys[2];
	b3Vec3 cap[B3_FLUID_MAX_FACES];
} b3ClipBuffers;

// Clips a polyhedron against a plane, keeping the submerged part behind it.
// Returns false if nothing is left.
static bool b3ClipPolyhedronByPlane( const b3ClipPolyhedron* in, b3Plane plane, b3ClipPolyhedron* out, b3ClipBuffers* buffers )
{
	out->pointCount = 0;
	out->faceCount = 0;

	b3Vec3* cap = buffers->cap;
	int capCount = 0;

	for ( int i = 0; i < in->faceCount; ++i )
	{
		const b3Vec3* points = in->points + in->faceStarts[i];
		int count = in->faceCounts[i];

		b3BeginFace( out );

		b3Vec3 a = points[count - 1];
		float da = b3PlaneSeparation( plane, a );

		for ( int j = 0; j < count; ++j )
		{
			b3Vec3 b = points[j];
			float db = b3PlaneSeparation( plane, b );

			if ( ( da <= 0.0f ) != ( db <= 0.0f ) )
			{
				// a-b crosses the plane
				float t = da / ( da - db );
				b3Vec3 p = b3Lerp( a, b, t );
				b3PushPoint( out, p );

				if ( capCount < B3_FLUID_MAX_FACES )
				{
					cap[capCount++] = p;
				}
			}

			if ( db <= 0.0f )
			{
				b3PushPoint( out, b );
			}

			a = b;
			da = db;
		}

		b3EndFace( out );
	}

	if ( capCount >= 3 )
	{
		// The cap lies on the plane.
		// Sort around its center so the winding is CCW about the plane normal
		b3Vec3 center = b3Vec3_zero;
		for ( int i = 0; i < capCount; ++i )
		{
			center = b3Add( center, cap[i] );
		}
		center = b3MulSV( 1.0f / capCount, center );

		b3Vec3 u = b3Perp( plane.normal );
		b3Vec3 v = b3Cross( plane.normal, u );

		// sort indices by angle
		b3CapSort sorted[B3_FLUID_MAX_FACES];
		for ( int i = 0; i < capCount; ++i )
		{
			b3Vec3 d = b3Sub( cap[i], center );
			sorted[i].angle = b3Atan2( b3Dot( d, v ), b3Dot( d, u ) );
			sorted[i].index = i;
		}

		qsort( sorted, capCount, sizeof( b3CapSort ), b3CompareCapAngles );

		b3BeginFace( out );
		for ( int i = 0; i < capCount; ++i )
		{
			b3PushPoint( out, cap[sorted[i].index] );
		}
		b3EndFace( out );
	}

	return out->faceCount >= 4;
}

static float b3TetrahedronVolume( b3Vec3* c, b3Vec3 v1, b3Vec3 v2, b3Vec3 v3 )
{
	float volume = b3Dot( v1, b3Cross( v2, v3 ) );
	*c = b3MulAdd( *c, volume, b3Add( v1, b3Add( v2, v3 ) ) );
	return volume;
}

// Computes the volume and centroid of a polyhedron
// https://people.eecs.berkeley.edu/~jfc/mirtich/massProps.html
static float b3ComputeVolume( const b3ClipPolyhedron* poly, b3Vec3* centroid )
{
	if ( poly->faceCount == 0 )
	{
		return 0.0f;
	}

	b3Vec3 p = poly->points[0];

	float volume = 0.0f;
	b3Vec3 c = b3Vec3_zero;

	// Compute the contribution of each face
	for ( int i = 0; i < poly->faceCount; ++i )
	{
		const b3Vec3* points = poly->points + poly->faceStarts[i];
		int count = poly->faceCounts[i];

		b3Vec3 v1 = b3Sub( points[0], p );
		for ( int j = 1; j + 1 < count; ++j )
		{
			b3Vec3 v2 = b3Sub( points[j], p );
			b3Vec3 v3 = b3Sub( points[j + 1], p );
			volume += b3TetrahedronVolume( &c, v1, v2, v3 );
		}
	}

	// Prevent slivers
	if ( volume <= FLT_EPSILON )
	{
		return 0.0f;
	}

	// Normalize
	*centroid = b3MulAdd( p, 0.25f / volume, c );

	return volume / 6.0f;
}

// Submerged part of a sphere behind a plane
// https://en.wikipedia.org/wiki/Spherical_cap
static float b3SphereCapFraction( b3Vec3 center, float radius, b3Plane plane, b3Vec3* offset )
{
	// height of the cap behind the plane
	float h = b3ClampFloat( radius - b3PlaneSeparation( plane, center ), 0.0f, 2.0f * radius );
	*offset = b3Vec3_zero;

	if ( h <= 0.0f )
	{
		return 0.0f;
	}

	if ( h >= 2.0f * radius )
	{
		return 1.0f;
	}

	// cap volume over the sphere volum
	float r3 = radius * radius * radius;
	float fraction = h * h * ( 3.0f * radius - h ) / ( 4.0f * r3 );

	// cap centroid distance from the center
	float d = 2.0f * radius - h;
	float distance = 3.0f * d * d / ( 4.0f * ( 3.0f * radius - h ) );

	*offset = b3MulSV( -distance, plane.normal );

	return fraction;
}

// approximation of the submerged part of a sphere against several planes
static float b3SphereSubmergedFraction( b3Vec3 center, float radius, const b3Plane* planes, int planeCount, b3Vec3* centroid )
{
	b3Vec3 offset;
	float fraction = b3SphereCapFraction( center, radius, planes[0], &offset );
	b3Vec3 total = offset;

	for ( int i = 1; i < planeCount && fraction > 0.0f; ++i )
	{
		float planeFraction = b3SphereCapFraction( center, radius, planes[i], &offset );
		if ( planeFraction < 1.0f )
		{
			fraction *= planeFraction;
			total = b3Add( total, offset );
		}
	}

	*centroid = b3Add( center, total );

	return fraction;
}

typedef struct b3Submerged
{
	float volume;
	b3Vec3 centroid;
	float length2;
} b3Submerged;

// Computes the submerged volume and center of buoyancy of a shape
static bool b3SubmergedVolume( const b3Shape* shape, b3Transform transform, const b3Plane* planes, int planeCount,
								b3ClipBuffers* buffers, b3Submerged* result )
{
	result->volume = 0.0f;
	result->centroid = b3Vec3_zero;
	result->length2 = 0.0f;

	switch ( shape->type )
	{
		case b3_sphereShape:
		{
			float radius = shape->sphere.radius;
			b3Vec3 center = b3TransformPoint( transform, shape->sphere.center );
			float fraction = b3SphereSubmergedFraction( center, radius, planes, planeCount, &result->centroid );
			result->volume = fraction * ( 4.0f / 3.0f ) * B3_PI * radius * radius * radius;
			result->length2 = 0.4f * radius * radius;
			return result->volume > 0.0f;
		}

		case b3_capsuleShape:
		{
			float radius = shape->capsule.radius;
			b3Vec3 p1 = b3TransformPoint( transform, shape->capsule.center1 );
			b3Vec3 p2 = b3TransformPoint( transform, shape->capsule.center2 );
			float length = b3Distance( p1, p2 );

			// treat the capsule as a row of spheres
			// TODO not happy with this approximation. need something better
			const int capsuleSamples = 5;

			float fraction = 0.0f;
			b3Vec3 sum = b3Vec3_zero;
			for ( int i = 0; i < capsuleSamples; ++i )
			{
				float t = (float)i / ( capsuleSamples - 1 );

				b3Vec3 centroid;
				float sampleFraction = b3SphereSubmergedFraction( b3Lerp( p1, p2, t ), radius, planes, planeCount, &centroid );

				fraction += sampleFraction;
				sum = b3MulAdd( sum, sampleFraction, centroid );
			}

			if ( fraction <= 0.0f )
			{
				return false;
			}

			float volume = B3_PI * radius * radius * ( length + ( 4.0f / 3.0f ) * radius );
			result->volume = volume * fraction / capsuleSamples;
			result->centroid = b3MulSV( 1.0f / fraction, sum );
			result->length2 = length * length / 12.0f + 0.4f * radius * radius;
			return true;
		}

		case b3_hullShape:
		{
			const b3HullData* hull = shape->hull;
			b3Vec3 extent = b3MulSV( 0.5f, b3Sub( hull->aabb.upperBound, hull->aabb.lowerBound ) );
			result->length2 = b3Dot( extent, extent ) / 3.0f;

			// compute the depth of each vertex, for a speedy reject
			const b3Vec3* points = b3GetHullPoints( hull );
			int vertexCount = hull->vertexCount;

			b3Plane clipPlanes[B3_FLUID_MAX_PLANES];
			int clipCount = 0;

			for ( int i = 0; i < planeCount; ++i )
			{
				// transform plane into the polyhedron frame.
				b3Plane plane = planes[i];
				b3Plane localPlane = {
					b3InvRotateVector( transform.q, plane.normal ),
					plane.offset - b3Dot( plane.normal, transform.p ),
				};

				float minDistance = FLT_MAX, maxDistance = -FLT_MAX;
				for ( int j = 0; j < vertexCount; ++j )
				{
					float distance = b3PlaneSeparation( localPlane, points[j] );
					minDistance = b3MinFloat( minDistance, distance );
					maxDistance = b3MaxFloat( maxDistance, distance );
				}

				if ( minDistance >= 0.0f )
				{
				    // Nothing is submerged
					return false;
				}

				if ( maxDistance > 0.0f )
				{
					clipPlanes[clipCount++] = plane;
				}
			}

			if ( clipCount == 0 )
			{
				// Fully submerged
				result->volume = hull->volume;
				result->centroid = b3TransformPoint( transform, hull->center );
				return hull->volume > 0.0f;
			}

			b3ClipPolyhedron* poly = buffers->polys + 0;
			b3ClipPolyhedron* clipped = buffers->polys + 1;
			if ( b3MakeHullPolyhedron( poly, hull, transform ) == false )
			{
				return false;
			}

			for ( int i = 0; i < clipCount; ++i )
			{
				if ( b3ClipPolyhedronByPlane( poly, clipPlanes[i], clipped, buffers ) == false )
				{
					return false;
				}

				b3ClipPolyhedron* temp = poly;
				poly = clipped;
				clipped = temp;
			}

			result->volume = b3ComputeVolume( poly, &result->centroid );
			// clipping can't add volume
			result->volume = b3MinFloat( result->volume, hull->volume );

			return result->volume > 0.0f;
		}

		default:
	    {
			return false;
		}
	}
}

// Gather the clip planes of a fluid sensor in the frame of the sensor body
static int b3GetFluidPlanes( const b3Shape* sensorShape, const b3Fluid* fluid, b3Plane* planes )
{
	int count = 0;
	planes[count++] = fluid->surface;

	if ( sensorShape->type == b3_hullShape )
	{
		const b3HullData* hull = sensorShape->hull;
		const b3Plane* hullPlanes = b3GetHullPlanes( hull );
		for ( int i = 0; i < hull->faceCount && count < B3_FLUID_MAX_PLANES; ++i )
		{
			b3Plane plane = hullPlanes[i];

			// Don't cut faces parallel to the surface
			if ( b3Dot( plane.normal, fluid->surface.normal ) > 0.999f )
			{
			    if ( plane.offset >= fluid->surface.offset - B3_LINEAR_SLOP )
				{
				    continue;
				}
			}

			planes[count++] = plane;
		}
	}

	return count;
}

static b3Vec3 b3GetFluidVelocity( b3World* world, b3Body* sensorBody, const b3Fluid* fluid, b3Pos base, b3Vec3 point )
{
	b3Vec3 velocity = fluid->currentVelocity;

	if ( sensorBody->setIndex == b3_awakeSet )
	{
		b3BodyState* state = b3GetBodyState( world, sensorBody );
		b3BodySim* sim = b3GetBodySim( world, sensorBody );

		b3Vec3 r = b3Sub( point, b3SubPos( sim->center, base ) );
		velocity = b3Add( velocity, b3Add( state->linearVelocity, b3Cross( state->angularVelocity, r ) ) );
	}

	return velocity;
}

// All the stuff one submerged shape does to its body
typedef struct b3FluidContribution
{
	int bodyId;
	int order;

	b3Vec3 leverArm;

	b3Vec3 buoyancyForce;
	b3Vec3 waterVelocity;

	float linearDrag;
	float angularDrag;
} b3FluidContribution;

static int b3CompareContributions( const void* a, const void* b )
{
	const b3FluidContribution* ca = (const b3FluidContribution*)a;
	const b3FluidContribution* cb = (const b3FluidContribution*)b;
	if ( ca->bodyId != cb->bodyId )
	{
		return ca->bodyId < cb->bodyId ? -1 : 1;
	}

	// break tie by order to keep sums deterministic
	return ( ca->order > cb->order ) - ( ca->order < cb->order );
}

static int b3GatherVisitors( b3World* world, b3Sensor* sensor, const b3Visitor* visitors, int visitorCount, b3Body* sensorBody,
							 b3Transform sensorTransform, b3Pos base, const b3Plane* planes, int planeCount, b3ClipBuffers* buffers,
							 b3FluidContribution* contributions, int order )
{
	const b3Fluid* fluid = &sensor->fluid;
	b3Vec3 gravity = world->gravity;
	int count = 0;

	for ( int i = 0; i < visitorCount; ++i )
	{
		b3Shape* shape = b3Array_Get( world->shapes, visitors[i].shapeId );
		if ( shape->id != visitors[i].shapeId || shape->generation != visitors[i].generation )
		{
			continue;
		}

		if ( shape->buoyancyScale == 0.0f && shape->fluidDragScale == 0.0f )
		{
			continue;
		}

		b3Body* body = b3Array_Get( world->bodies, shape->bodyId );
		if ( body->setIndex != b3_awakeSet || body->type != b3_dynamicBody )
		{
			continue;
		}

		b3BodySim* sim = b3GetBodySim( world, body );
		if ( sim->invMass == 0.0f )
		{
			continue;
		}

		// get visitor shape in the fluid frame
		b3Transform bodyTransform = b3ToRelativeTransform( sim->transform, base );
		b3Transform relativeTransform = b3InvMulTransforms( sensorTransform, bodyTransform );

		b3Submerged submerged;
		if ( b3SubmergedVolume( shape, relativeTransform, planes, planeCount, buffers, &submerged ) == false )
		{
			continue;
		}

		b3Vec3 point = b3TransformPoint( sensorTransform, submerged.centroid );
		float displacedMass = fluid->density * submerged.volume;

		float partialMass = displacedMass * shape->fluidDragScale;

		b3FluidContribution* c = contributions + count;
		c->bodyId = shape->bodyId;
		c->order = order + count;
		c->leverArm = b3Sub( point, b3SubPos( sim->center, base ) );

		// buoyancy = -density * volume * gravity
		c->buoyancyForce = b3MulSV( -displacedMass * shape->buoyancyScale, gravity );
		c->waterVelocity = b3GetFluidVelocity( world, sensorBody, fluid, base, point );

		// quadratic term
		b3BodyState* state = b3GetBodyState( world, body );
		b3Vec3 vc = b3Add( state->linearVelocity, b3Cross( state->angularVelocity, c->leverArm ) );
		float speed = b3Length( b3Sub( vc, c->waterVelocity ) );

		c->linearDrag = ( fluid->linearDrag + fluid->quadraticDrag * speed ) * partialMass;
		c->angularDrag = fluid->angularDrag * partialMass * submerged.length2;

		count += 1;
	}

	return count;
}

static int b3GatherFluid( b3World* world, b3Sensor* sensor, b3ClipBuffers* buffers, b3FluidContribution* contributions, int order )
{
	const b3Fluid* fluid = &sensor->fluid;
	b3Shape* sensorShape = b3Array_Get( world->shapes, sensor->shapeId );
	b3Body* sensorBody = b3Array_Get( world->bodies, sensorShape->bodyId );

	if ( sensorBody->setIndex == b3_disabledSet )
	{
		return 0;
	}

	// Work relative to the sensor body for precision
	b3WorldTransform sensorWorldTransform = b3GetBodyTransformQuick( world, sensorBody );
	b3Pos base = sensorWorldTransform.p;
	b3Transform sensorTransform = b3ToRelativeTransform( sensorWorldTransform, base );

	b3Plane planes[B3_FLUID_MAX_PLANES];
	int planeCount = b3GetFluidPlanes( sensorShape, fluid, planes );

	int count = b3GatherVisitors( world, sensor, sensor->overlaps2.data, sensor->overlaps2.count, sensorBody, sensorTransform, base,
								  planes, planeCount, buffers, contributions, order );
	count += b3GatherVisitors( world, sensor, sensor->fluidVisitors.data, sensor->fluidVisitors.count, sensorBody, sensorTransform,
							   base, planes, planeCount, buffers, contributions + count, order + count );

	return count;
}

// r x (w x r) = (|r|^2 I - r r^T) w
static b3Matrix3 b3OffsetDampingMatrix( b3Vec3 r, float k )
{
	float rr = b3Dot( r, r );
	b3Matrix3 m;
	m.cx = (b3Vec3){ k * ( rr - r.x * r.x ), -k * r.y * r.x, -k * r.z * r.x };
	m.cy = (b3Vec3){ -k * r.x * r.y, k * ( rr - r.y * r.y ), -k * r.z * r.y };
	m.cz = (b3Vec3){ -k * r.x * r.z, -k * r.y * r.z, k * ( rr - r.z * r.z ) };
	return m;
}

// Solve for implicit drag
static void b3ApplyFluidToBody( b3World* world, const b3FluidContribution* contributions, int count, float timeStep )
{
	b3Body* body = b3Array_Get( world->bodies, contributions[0].bodyId );
	b3BodySim* sim = b3GetBodySim( world, body );
	b3BodyState* state = b3GetBodyState( world, body );
	float h = timeStep;

	// linear: v' = v - h invMass sum k_i (v' - u_i)
	// angular: w' = w - h invI D w'

	float linearDrag = 0.0f;
	b3Vec3 dragTarget = b3Vec3_zero;

	b3Matrix3 damping = b3Mat3_zero;

	for ( int i = 0; i < count; ++i )
	{
		const b3FluidContribution* c = contributions + i;

		linearDrag += c->linearDrag;

		dragTarget = b3MulAdd( dragTarget, c->linearDrag, c->waterVelocity );

		damping = b3AddMM( damping, b3OffsetDampingMatrix( c->leverArm, c->linearDrag ) );
		damping.cx.x += c->angularDrag;
		damping.cy.y += c->angularDrag;
		damping.cz.z += c->angularDrag;
	}

	float hm = h * sim->invMass;
	b3Vec3 v = b3MulSV( 1.0f / ( 1.0f + hm * linearDrag ), b3MulAdd( state->linearVelocity, hm, dragTarget ) );

	b3Matrix3 system = b3AddMM( b3Mat3_identity, b3MulSM( h, b3MulMM( sim->invInertiaWorld, damping ) ) );
	b3Vec3 w = b3MulMV( b3InvertMatrix( system ), state->angularVelocity );

	b3Vec3 force = b3Vec3_zero;
	b3Vec3 torque = b3Neg( b3MulMV( damping, w ) );
	for ( int i = 0; i < count; ++i )
	{
		const b3FluidContribution* c = contributions + i;

		b3Vec3 totalForce = b3MulSub( c->buoyancyForce, c->linearDrag, b3Sub( v, c->waterVelocity ) );
		b3Vec3 spin = b3MulSV( -c->linearDrag, b3Cross( w, c->leverArm ) );

		force = b3Add( force, b3Add( totalForce, spin ) );
		torque = b3Add( torque, b3Cross( c->leverArm, totalForce ) );
	}

	sim->force = b3Add( sim->force, force );
	sim->torque = b3Add( sim->torque, torque );
}

static bool b3FluidBodyMoving( b3World* world, b3Body* body )
{
	if ( body->setIndex != b3_awakeSet )
	{
		return false;
	}

	b3BodyState* state = b3GetBodyState( world, body );
	b3BodySim* sim = b3GetBodySim( world, body );
	float speed = b3Length( state->linearVelocity ) + b3Length( state->angularVelocity ) * b3Length( sim->maxExtent );
	return speed > body->sleepThreshold;
}

static void b3WakeVisitors( b3World* world, const b3Visitor* visitors, int count )
{
	for ( int j = 0; j < count; ++j )
	{
		const b3Visitor* visitor = visitors + j;
		b3Shape* shape = b3Array_Get( world->shapes, visitor->shapeId );
		if ( shape->id != visitor->shapeId || shape->generation != visitor->generation )
		{
			continue;
		}

		b3Body* body = b3Array_Get( world->bodies, shape->bodyId );
		if ( body->type == b3_dynamicBody )
		{
			b3WakeBody( world, body );
		}
	}
}

void b3WakeFluidVisitors( b3World* world )
{
	int sensorCount = world->sensors.count;
	for ( int i = 0; i < sensorCount; ++i )
	{
		b3Sensor* sensor = world->sensors.data + i;
		bool wake = sensor->fluid.wakeVisitors;
		sensor->fluid.wakeVisitors = false;

		if ( sensor->fluid.enabled == false && wake == false )
		{
			continue;
		}

		b3Shape* sensorShape = b3Array_Get( world->shapes, sensor->shapeId );
		b3Body* sensorBody = b3Array_Get( world->bodies, sensorShape->bodyId );

		if ( wake == false && b3FluidBodyMoving( world, sensorBody ) == false )
		{
			continue;
		}

		b3WakeVisitors( world, sensor->overlaps2.data, sensor->overlaps2.count );
		b3WakeVisitors( world, sensor->fluidVisitors.data, sensor->fluidVisitors.count );
	}
}

void b3ApplyFluidForces( b3World* world, float timeStep )
{
	if ( timeStep <= 0.0f )
	{
		return;
	}

	b3TracyCZoneNC( fluids, "Fluids", b3_colorDodgerBlue, true );

	int maxCount = 0;
	int sensorCount = world->sensors.count;
	for ( int i = 0; i < sensorCount; ++i )
	{
		b3Sensor* sensor = world->sensors.data + i;
		if ( sensor->fluid.enabled )
		{
			maxCount += sensor->overlaps2.count + sensor->fluidVisitors.count;
		}
	}

	if ( maxCount == 0 )
	{
		b3TracyCZoneEnd( fluids );
		return;
	}

	b3FluidContribution* contributions = b3StackAlloc( &world->stack, maxCount * (int)sizeof( b3FluidContribution ), "fluid contributions" );
	b3ClipBuffers* buffers = b3StackAlloc( &world->stack, (int)sizeof( b3ClipBuffers ), "fluid clip" );

	int count = 0;
	for ( int i = 0; i < sensorCount; ++i )
	{
		b3Sensor* sensor = world->sensors.data + i;
		if ( sensor->fluid.enabled )
		{
			count += b3GatherFluid( world, sensor, buffers, contributions + count, count );
		}
	}

	qsort( contributions, count, sizeof( b3FluidContribution ), b3CompareContributions );

	int start = 0;
	while ( start < count )
	{
		int end = start + 1;
		while ( end < count && contributions[end].bodyId == contributions[start].bodyId )
		{
			end += 1;
		}

		b3ApplyFluidToBody( world, contributions + start, end - start, timeStep );
		start = end;
	}

	b3StackFree( &world->stack, buffers );
	b3StackFree( &world->stack, contributions );

	b3TracyCZoneEnd( fluids );
}

static b3Sensor* b3GetFluidSensor( b3World* world, b3Shape* shape )
{
	if ( shape->sensorIndex == B3_NULL_INDEX )
	{
		return NULL;
	}

	return b3Array_Get( world->sensors, shape->sensorIndex );
}

static b3Sensor* b3GetFluid( b3World* world, b3ShapeId shapeId )
{
	b3Shape* shape = b3GetShape( world, shapeId );
	b3Sensor* sensor = b3GetFluidSensor( world, shape );
	B3_ASSERT( sensor != NULL && sensor->fluid.enabled );
	return sensor != NULL && sensor->fluid.enabled ? sensor : NULL;
}

void b3Shape_EnableFluid( b3ShapeId shapeId, const b3FluidDef* def )
{
	B3_CHECK_DEF( def );
	B3_ASSERT( b3IsValidPlane( def->surface ) );
	B3_ASSERT( b3IsValidFloat( def->density ) && def->density >= 0.0f );
	B3_ASSERT( b3IsValidFloat( def->linearDrag ) && def->linearDrag >= 0.0f );
	B3_ASSERT( b3IsValidFloat( def->quadraticDrag ) && def->quadraticDrag >= 0.0f );
	B3_ASSERT( b3IsValidFloat( def->angularDrag ) && def->angularDrag >= 0.0f );
	B3_ASSERT( b3IsValidVec3( def->currentVelocity ) );

	b3World* world = b3GetUnlockedWorld( shapeId.world0 );
	if ( world == NULL )
	{
		return;
	}

	B3_REC( world, ShapeEnableFluid, shapeId, *def );

	b3Shape* shape = b3GetShape( world, shapeId );
	b3Sensor* sensor = b3GetFluidSensor( world, shape );
	B3_ASSERT( sensor != NULL );
	if ( sensor == NULL )
	{
		return;
	}

	sensor->fluid = (b3Fluid){
		.surface = def->surface,
		.currentVelocity = def->currentVelocity,
		.density = def->density,
		.linearDrag = def->linearDrag,
		.quadraticDrag = def->quadraticDrag,
		.angularDrag = def->angularDrag,
		.enabled = true,
		.wakeVisitors = true,
	};
}

void b3Shape_DisableFluid( b3ShapeId shapeId )
{
	b3World* world = b3GetUnlockedWorld( shapeId.world0 );
	if ( world == NULL )
	{
		return;
	}

	B3_REC( world, ShapeDisableFluid, shapeId );

	b3Sensor* sensor = b3GetFluid( world, shapeId );
	if ( sensor == NULL )
	{
		return;
	}

	sensor->fluid.enabled = false;
	sensor->fluid.wakeVisitors = true;
}

bool b3Shape_IsFluid( b3ShapeId shapeId )
{
	b3World* world = b3GetWorld( shapeId.world0 );
	b3Shape* shape = b3GetShape( world, shapeId );
	b3Sensor* sensor = b3GetFluidSensor( world, shape );
	return sensor != NULL && sensor->fluid.enabled;
}

b3FluidDef b3Shape_GetFluid( b3ShapeId shapeId )
{
	b3World* world = b3GetWorld( shapeId.world0 );
	b3Sensor* sensor = b3GetFluid( world, shapeId );

	b3FluidDef def = b3DefaultFluidDef();
	if ( sensor != NULL )
	{
		def.surface = sensor->fluid.surface;
		def.currentVelocity = sensor->fluid.currentVelocity;
		def.density = sensor->fluid.density;
		def.linearDrag = sensor->fluid.linearDrag;
		def.quadraticDrag = sensor->fluid.quadraticDrag;
		def.angularDrag = sensor->fluid.angularDrag;
	}

	return def;
}

void b3Shape_SetFluidSurface( b3ShapeId shapeId, b3Plane surface )
{
	B3_ASSERT( b3IsValidPlane( surface ) );

	b3World* world = b3GetUnlockedWorld( shapeId.world0 );
	if ( world == NULL )
	{
		return;
	}

	B3_REC( world, ShapeSetFluidSurface, shapeId, surface.normal, surface.offset );

	b3Sensor* sensor = b3GetFluid( world, shapeId );
	if ( sensor != NULL )
	{
		sensor->fluid.surface = surface;
		sensor->fluid.wakeVisitors = true;
	}
}

void b3Shape_SetFluidCurrent( b3ShapeId shapeId, b3Vec3 currentVelocity )
{
	B3_ASSERT( b3IsValidVec3( currentVelocity ) );

	b3World* world = b3GetUnlockedWorld( shapeId.world0 );
	if ( world == NULL )
	{
		return;
	}

	B3_REC( world, ShapeSetFluidCurrent, shapeId, currentVelocity );

	b3Sensor* sensor = b3GetFluid( world, shapeId );
	if ( sensor != NULL )
	{
		sensor->fluid.currentVelocity = currentVelocity;
		sensor->fluid.wakeVisitors = true;
	}
}

void b3Shape_SetBuoyancyScale( b3ShapeId shapeId, float buoyancyScale )
{
	B3_ASSERT( b3IsValidFloat( buoyancyScale ) && buoyancyScale >= 0.0f );

	b3World* world = b3GetUnlockedWorld( shapeId.world0 );
	if ( world == NULL )
	{
		return;
	}

	B3_REC( world, ShapeSetBuoyancyScale, shapeId, buoyancyScale );

	b3Shape* shape = b3GetShape( world, shapeId );
	shape->buoyancyScale = buoyancyScale;
}

float b3Shape_GetBuoyancyScale( b3ShapeId shapeId )
{
	b3World* world = b3GetWorld( shapeId.world0 );
	b3Shape* shape = b3GetShape( world, shapeId );
	return shape->buoyancyScale;
}

void b3Shape_SetFluidDragScale( b3ShapeId shapeId, float fluidDragScale )
{
	B3_ASSERT( b3IsValidFloat( fluidDragScale ) && fluidDragScale >= 0.0f );

	b3World* world = b3GetUnlockedWorld( shapeId.world0 );
	if ( world == NULL )
	{
		return;
	}

	B3_REC( world, ShapeSetFluidDragScale, shapeId, fluidDragScale );

	b3Shape* shape = b3GetShape( world, shapeId );
	shape->fluidDragScale = fluidDragScale;
}

float b3Shape_GetFluidDragScale( b3ShapeId shapeId )
{
	b3World* world = b3GetWorld( shapeId.world0 );
	b3Shape* shape = b3GetShape( world, shapeId );
	return shape->fluidDragScale;
}

bool b3Shape_GetSubmergedVolume( b3ShapeId shapeId, b3ShapeId fluidShapeId, float* volume, b3Pos* centroid )
{
	*volume = 0.0f;

	b3World* world = b3GetUnlockedWorld( shapeId.world0 );
	if ( world == NULL )
	{
		return false;
	}

	b3Shape* shape = b3GetShape( world, shapeId );
	b3Shape* sensorShape = b3GetShape( world, fluidShapeId );
	b3Sensor* sensor = b3GetFluid( world, fluidShapeId );
	if ( sensor == NULL )
	{
		return false;
	}

	b3WorldTransform sensorWorldTransform = b3GetBodyTransform( world, sensorShape->bodyId );
	b3Pos base = sensorWorldTransform.p;
	b3Transform sensorTransform = b3ToRelativeTransform( sensorWorldTransform, base );
	b3Transform bodyTransform = b3ToRelativeTransform( b3GetBodyTransform( world, shape->bodyId ), base );
	b3Transform relativeTransform = b3InvMulTransforms( sensorTransform, bodyTransform );

	b3Plane planes[B3_FLUID_MAX_PLANES];
	int planeCount = b3GetFluidPlanes( sensorShape, &sensor->fluid, planes );

	b3ClipBuffers* buffers = b3StackAlloc( &world->stack, (int)sizeof( b3ClipBuffers ), "fluid clip" );
	b3Submerged submerged;
	bool result = b3SubmergedVolume( shape, relativeTransform, planes, planeCount, buffers, &submerged );
	b3StackFree( &world->stack, buffers );

	if ( result == false )
	{
		return false;
	}

	*volume = submerged.volume;
	if ( centroid != NULL )
	{
		*centroid = b3OffsetPos( base, b3TransformPoint( sensorTransform, submerged.centroid ) );
	}

	return true;
}
