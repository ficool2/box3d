// SPDX-FileCopyrightText: 2026 ficool2
// SPDX-License-Identifier: MIT

#include "test_macros.h"

#include "box3d/box3d.h"
#include "box3d/collision.h"
#include "box3d/math_functions.h"

#include <float.h>
#include <math.h>

// A pool: 20 x 20 x 10 sensor box with the surface at z = 0
static b3ShapeId CreatePool( b3WorldId worldId, b3BoxHull* poolHull )
{
	*poolHull = b3MakeOffsetBoxHull( 10.0f, 10.0f, 5.0f, (b3Vec3){ 0.0f, 0.0f, -5.0f } );

	b3BodyDef bodyDef = b3DefaultBodyDef();
	b3BodyId bodyId = b3CreateBody( worldId, &bodyDef );

	b3ShapeDef shapeDef = b3DefaultShapeDef();
	shapeDef.isSensor = true;
	shapeDef.enableSensorEvents = true;
	b3ShapeId shapeId = b3CreateHullShape( bodyId, &shapeDef, &poolHull->base );

	b3FluidDef fluidDef = b3DefaultFluidDef();
	fluidDef.surface = (b3Plane){ { 0.0f, 0.0f, 1.0f }, 0.0f };
	fluidDef.density = 1000.0f;
	b3Shape_EnableFluid( shapeId, &fluidDef );
	return shapeId;
}

static b3ShapeId CreateBox( b3WorldId worldId, const b3BoxHull* hull, b3Vec3 position, b3Quat rotation, float density )
{
	b3BodyDef bodyDef = b3DefaultBodyDef();
	bodyDef.type = b3_dynamicBody;
	bodyDef.position = b3ToPos( position );
	bodyDef.rotation = rotation;
	b3BodyId bodyId = b3CreateBody( worldId, &bodyDef );

	b3ShapeDef shapeDef = b3DefaultShapeDef();
	shapeDef.density = density;
	shapeDef.enableSensorEvents = true;
	return b3CreateHullShape( bodyId, &shapeDef, &hull->base );
}

static int SubmergedVolumeTest( void )
{
	b3WorldDef worldDef = b3DefaultWorldDef();
	b3WorldId worldId = b3CreateWorld( &worldDef );

	b3BoxHull poolHull;
	b3ShapeId poolId = CreatePool( worldId, &poolHull );
	ENSURE( b3Shape_IsFluid( poolId ) );

	b3BoxHull cube = b3MakeCubeHull( 0.5f );
	float volume;
	b3Pos centroid;

	// Half under the surface
	b3ShapeId boxId = CreateBox( worldId, &cube, b3Vec3_zero, b3Quat_identity, 1000.0f );
	ENSURE( b3Shape_GetSubmergedVolume( boxId, poolId, &volume, &centroid ) );
	ENSURE_SMALL( volume - 0.5f, 1e-5f );
	ENSURE_SMALL( centroid.z + 0.25f, 1e-5f );
	ENSURE_SMALL( centroid.x, 1e-5f );

	// Rotated 45 degrees about x
	b3BodyId bodyId = b3Shape_GetBody( boxId );
	b3Quat q = b3MakeQuatFromAxisAngle( (b3Vec3){ 1.0f, 0.0f, 0.0f }, 0.25f * B3_PI );
	b3Body_SetTransform( bodyId, b3Pos_zero, q );
	ENSURE( b3Shape_GetSubmergedVolume( boxId, poolId, &volume, &centroid ) );
	ENSURE_SMALL( volume - 0.5f, 1e-5f );
	ENSURE_SMALL( centroid.z + sqrtf( 0.5f ) / 3.0f, 1e-3f );

	// Fully submerged
	b3Body_SetTransform( bodyId, b3ToPos( (b3Vec3){ 0.0f, 0.0f, -2.0f } ), b3Quat_identity );
	ENSURE( b3Shape_GetSubmergedVolume( boxId, poolId, &volume, &centroid ) );
	ENSURE_SMALL( volume - 1.0f, 1e-5f );

	// Straddling the side wall of the pool at x = 10
	b3Body_SetTransform( bodyId, b3ToPos( (b3Vec3){ 10.0f, 0.0f, 0.0f } ), b3Quat_identity );
	ENSURE( b3Shape_GetSubmergedVolume( boxId, poolId, &volume, &centroid ) );
	ENSURE_SMALL( volume - 0.25f, 1e-5f );
	ENSURE_SMALL( centroid.x - 9.75f, 1e-4f );

	// Above the water
	b3Body_SetTransform( bodyId, b3ToPos( (b3Vec3){ 0.0f, 0.0f, 2.0f } ), b3Quat_identity );
	ENSURE( b3Shape_GetSubmergedVolume( boxId, poolId, &volume, &centroid ) == false );

	// Sphere half under
	{
		b3BodyDef bodyDef = b3DefaultBodyDef();
		bodyDef.type = b3_dynamicBody;
		b3BodyId sphereBodyId = b3CreateBody( worldId, &bodyDef );
		b3ShapeDef shapeDef = b3DefaultShapeDef();
		b3Sphere sphere = { b3Vec3_zero, 1.0f };
		b3ShapeId sphereId = b3CreateSphereShape( sphereBodyId, &shapeDef, &sphere );
		ENSURE( b3Shape_GetSubmergedVolume( sphereId, poolId, &volume, &centroid ) );
		ENSURE_SMALL( volume - ( 2.0f / 3.0f ) * B3_PI, 1e-4f );
		ENSURE_SMALL( centroid.z + 3.0f / 8.0f, 1e-5f );
	}

	b3DestroyWorld( worldId );
	return 0;
}

static float SimulateBox( float density, float* outSpeed )
{
	b3WorldDef worldDef = b3DefaultWorldDef();
	worldDef.gravity = (b3Vec3){ 0.0f, 0.0f, -10.0f };
	b3WorldId worldId = b3CreateWorld( &worldDef );

	b3BoxHull poolHull;
	CreatePool( worldId, &poolHull );

	b3BoxHull cube = b3MakeCubeHull( 0.5f );
	b3ShapeId boxId = CreateBox( worldId, &cube, (b3Vec3){ 0.0f, 0.0f, 1.0f }, b3Quat_identity, density );
	b3BodyId bodyId = b3Shape_GetBody( boxId );

	for ( int i = 0; i < 600; ++i )
	{
		b3World_Step( worldId, 1.0f / 60.0f, 4 );
	}

	b3Pos p = b3Body_GetPosition( bodyId );
	*outSpeed = b3Length( b3Body_GetLinearVelocity( bodyId ) );
	b3DestroyWorld( worldId );
	return (float)p.z;
}

static int BuoyancyTest( void )
{
	float speed;

	// Half the density of water floats half submerged
	float z = SimulateBox( 500.0f, &speed );
	ENSURE_SMALL( z, 0.02f );
	ENSURE( speed < 0.05f );

	// A quarter of the density floats higher
	z = SimulateBox( 250.0f, &speed );
	ENSURE_SMALL( z - 0.25f, 0.02f );

	// Denser than water sinks to the bottom of the pool
	z = SimulateBox( 2000.0f, &speed );
	ENSURE( z < -10.0f );

	return 0;
}

static int FluidDragTest( void )
{
	b3WorldDef worldDef = b3DefaultWorldDef();
	worldDef.gravity = b3Vec3_zero;
	b3WorldId worldId = b3CreateWorld( &worldDef );

	b3BoxHull poolHull;
	b3ShapeId poolId = CreatePool( worldId, &poolHull );

	b3BoxHull cube = b3MakeCubeHull( 0.5f );
	b3ShapeId boxId = CreateBox( worldId, &cube, (b3Vec3){ 0.0f, 0.0f, -3.0f }, b3Quat_identity, 1000.0f );
	b3BodyId bodyId = b3Shape_GetBody( boxId );
	b3Body_SetLinearVelocity( bodyId, (b3Vec3){ 1.0f, 0.0f, 0.0f } );
	b3Body_SetAngularVelocity( bodyId, (b3Vec3){ 0.0f, 0.0f, 2.0f } );

	// one step to find the overlap
	b3World_Step( worldId, 1.0f / 60.0f, 4 );

	for ( int i = 0; i < 60; ++i )
	{
		b3World_Step( worldId, 1.0f / 60.0f, 4 );
	}

	// drag rate 1/s
	b3Vec3 v = b3Body_GetLinearVelocity( bodyId );
	ENSURE( v.x > 0.25f && v.x < 0.45f );

	b3Vec3 w = b3Body_GetAngularVelocity( bodyId );
	ENSURE( w.z > 0.0f && w.z < 2.0f );

	// Current drags the body along
	b3FluidDef fluidDef = b3Shape_GetFluid( poolId );
	fluidDef.currentVelocity = (b3Vec3){ 0.0f, 2.0f, 0.0f };
	b3Shape_EnableFluid( poolId, &fluidDef );
	for ( int i = 0; i < 300; ++i )
	{
		b3World_Step( worldId, 1.0f / 60.0f, 4 );
	}
	v = b3Body_GetLinearVelocity( bodyId );
	ENSURE_SMALL( v.y - 2.0f, 0.05f );

	// No drag when disabled
	b3Shape_SetFluidDragScale( boxId, 0.0f );
	b3Body_SetLinearVelocity( bodyId, (b3Vec3){ 1.0f, 0.0f, 0.0f } );
	b3World_Step( worldId, 1.0f / 60.0f, 4 );
	v = b3Body_GetLinearVelocity( bodyId );
	ENSURE_SMALL( v.x - 1.0f, 1e-5f );

	b3Shape_DisableFluid( poolId );
	ENSURE( b3Shape_IsFluid( poolId ) == false );

	b3DestroyWorld( worldId );
	return 0;
}

// A very light body made of many small hulls
static int MultiShapeDragTest( void )
{
	b3WorldDef worldDef = b3DefaultWorldDef();
	worldDef.gravity = b3Vec3_zero;
	b3WorldId worldId = b3CreateWorld( &worldDef );

	b3BoxHull poolHull;
	b3ShapeId poolId = CreatePool( worldId, &poolHull );

	b3BodyDef bodyDef = b3DefaultBodyDef();
	bodyDef.type = b3_dynamicBody;
	bodyDef.position = b3ToPos( (b3Vec3){ 0.0f, 0.0f, -3.0f } );
	b3BodyId bodyId = b3CreateBody( worldId, &bodyDef );

	b3ShapeDef shapeDef = b3DefaultShapeDef();
	shapeDef.density = 10.0f;
	shapeDef.enableSensorEvents = true;
	b3BoxHull pieces[8];
	for ( int i = 0; i < 8; ++i )
	{
		b3Vec3 offset = { ( i & 1 ) ? 0.25f : -0.25f, ( i & 2 ) ? 0.25f : -0.25f, ( i & 4 ) ? 0.25f : -0.25f };
		pieces[i] = b3MakeOffsetBoxHull( 0.25f, 0.25f, 0.25f, offset );
		b3CreateHullShape( bodyId, &shapeDef, &pieces[i].base );
	}

	b3FluidDef fluidDef = b3Shape_GetFluid( poolId );
	fluidDef.linearDrag = 4.0f;
	b3Shape_EnableFluid( poolId, &fluidDef );

	b3World_Step( worldId, 1.0f / 60.0f, 4 );
	b3Body_SetLinearVelocity( bodyId, (b3Vec3){ 1.0f, 0.0f, 0.0f } );

	for ( int i = 0; i < 10; ++i )
	{
		b3World_Step( worldId, 1.0f / 60.0f, 4 );
		b3Vec3 v = b3Body_GetLinearVelocity( bodyId );
		ENSURE( v.x >= -1e-3f && v.x <= 1.0f );
	}

	b3DestroyWorld( worldId );
	return 0;
}

// Changing the fluid wakes the bodies sleeping in it
static int FluidWakeTest( void )
{
	b3WorldDef worldDef = b3DefaultWorldDef();
	worldDef.gravity = (b3Vec3){ 0.0f, 0.0f, -10.0f };
	b3WorldId worldId = b3CreateWorld( &worldDef );

	b3BoxHull poolHull;
	b3ShapeId poolId = CreatePool( worldId, &poolHull );

	b3BoxHull cube = b3MakeCubeHull( 0.5f );
	b3ShapeId boxId = CreateBox( worldId, &cube, (b3Vec3){ 0.0f, 0.0f, 0.5f }, b3Quat_identity, 500.0f );
	b3BodyId bodyId = b3Shape_GetBody( boxId );

	for ( int i = 0; i < 600 && b3Body_IsAwake( bodyId ); ++i )
	{
		b3World_Step( worldId, 1.0f / 60.0f, 4 );
	}
	ENSURE( b3Body_IsAwake( bodyId ) == false );

	// lower the water
	b3Shape_SetFluidSurface( poolId, (b3Plane){ { 0.0f, 0.0f, 1.0f }, -0.25f } );
	b3World_Step( worldId, 1.0f / 60.0f, 4 );

	// the box wakes
	ENSURE( b3Body_IsAwake( bodyId ) );

	// and sinks
	for ( int i = 0; i < 600; ++i )
	{
		b3World_Step( worldId, 1.0f / 60.0f, 4 );
	}
	ENSURE_SMALL( (float)b3Body_GetPosition( bodyId ).z + 0.25f, 0.02f );

	b3DestroyWorld( worldId );
	return 0;
}

// Fluid finds the box with no sensor events anywhere
static int NoEventsTest( void )
{
	b3WorldDef worldDef = b3DefaultWorldDef();
	worldDef.gravity = (b3Vec3){ 0.0f, 0.0f, -10.0f };
	b3WorldId worldId = b3CreateWorld( &worldDef );

	b3BoxHull poolHull = b3MakeOffsetBoxHull( 10.0f, 10.0f, 5.0f, (b3Vec3){ 0.0f, 0.0f, -5.0f } );
	b3BodyDef poolBodyDef = b3DefaultBodyDef();
	b3BodyId poolBodyId = b3CreateBody( worldId, &poolBodyDef );
	b3ShapeDef poolDef = b3DefaultShapeDef();
	poolDef.isSensor = true;
	b3ShapeId poolId = b3CreateHullShape( poolBodyId, &poolDef, &poolHull.base );
	b3FluidDef fluidDef = b3DefaultFluidDef();
	fluidDef.density = 1000.0f;
	b3Shape_EnableFluid( poolId, &fluidDef );

	b3BoxHull cube = b3MakeCubeHull( 0.5f );
	b3BodyDef bodyDef = b3DefaultBodyDef();
	bodyDef.type = b3_dynamicBody;
	bodyDef.position = b3ToPos( (b3Vec3){ 0.0f, 0.0f, 1.0f } );
	b3BodyId bodyId = b3CreateBody( worldId, &bodyDef );
	b3ShapeDef shapeDef = b3DefaultShapeDef();
	shapeDef.density = 500.0f;
	b3CreateHullShape( bodyId, &shapeDef, &cube.base );

	for ( int i = 0; i < 600; ++i )
	{
		b3World_Step( worldId, 1.0f / 60.0f, 4 );
	}

	ENSURE_SMALL( (float)b3Body_GetPosition( bodyId ).z, 0.02f );
	ENSURE( b3World_GetSensorEvents( worldId ).beginCount == 0 );

	b3DestroyWorld( worldId );
	return 0;
}

// A long, very light rod spinning fast
static int LongRodTest( void )
{
	b3WorldDef worldDef = b3DefaultWorldDef();
	worldDef.gravity = b3Vec3_zero;
	b3WorldId worldId = b3CreateWorld( &worldDef );

	b3BoxHull poolHull;
	CreatePool( worldId, &poolHull );

	b3BoxHull rod = b3MakeBoxHull( 4.0f, 0.05f, 0.05f );
	b3ShapeId rodId = CreateBox( worldId, &rod, (b3Vec3){ 0.0f, 0.0f, -3.0f }, b3Quat_identity, 1.0f );
	b3BodyId bodyId = b3Shape_GetBody( rodId );
	b3World_Step( worldId, 1.0f / 60.0f, 4 );
	b3Body_SetAngularVelocity( bodyId, (b3Vec3){ 0.0f, 0.0f, 10.0f } );

	float last = 10.0f;
	for ( int i = 0; i < 60; ++i )
	{
		b3World_Step( worldId, 1.0f / 60.0f, 4 );
		b3Vec3 w = b3Body_GetAngularVelocity( bodyId );
		ENSURE( w.z >= -1e-3f && w.z <= last + 1e-3f );
		last = w.z;
	}

	b3DestroyWorld( worldId );
	return 0;
}

#if defined( BOX3D_DOUBLE_PRECISION )
// A pool far from the origin
static int FarFluidTest( void )
{
	b3WorldDef worldDef = b3DefaultWorldDef();
	worldDef.gravity = (b3Vec3){ 0.0f, 0.0f, -10.0f };
	b3WorldId worldId = b3CreateWorld( &worldDef );

	b3Pos far = { 1.0e7, -2.0e7, 5.0e6 };

	b3BoxHull poolHull = b3MakeOffsetBoxHull( 10.0f, 10.0f, 5.0f, (b3Vec3){ 0.0f, 0.0f, -5.0f } );
	b3BodyDef poolBodyDef = b3DefaultBodyDef();
	poolBodyDef.position = far;
	b3BodyId poolBodyId = b3CreateBody( worldId, &poolBodyDef );
	b3ShapeDef poolDef = b3DefaultShapeDef();
	poolDef.isSensor = true;
	b3ShapeId poolId = b3CreateHullShape( poolBodyId, &poolDef, &poolHull.base );
	b3FluidDef fluidDef = b3DefaultFluidDef();
	fluidDef.density = 1000.0f;
	b3Shape_EnableFluid( poolId, &fluidDef );

	b3BoxHull cube = b3MakeCubeHull( 0.5f );
	b3BodyDef bodyDef = b3DefaultBodyDef();
	bodyDef.type = b3_dynamicBody;
	bodyDef.position = b3OffsetPos( far, (b3Vec3){ 0.0f, 0.0f, 1.0f } );
	b3BodyId bodyId = b3CreateBody( worldId, &bodyDef );
	b3ShapeDef shapeDef = b3DefaultShapeDef();
	shapeDef.density = 500.0f;
	b3ShapeId boxId = b3CreateHullShape( bodyId, &shapeDef, &cube.base );

	float volume;
	b3Pos centroid;
	b3Body_SetTransform( bodyId, far, b3Quat_identity );
	ENSURE( b3Shape_GetSubmergedVolume( boxId, poolId, &volume, &centroid ) );
	ENSURE_SMALL( volume - 0.5f, 1e-5f );
	ENSURE_SMALL( (float)( centroid.z - far.z ) + 0.25f, 1e-5f );
	b3Body_SetTransform( bodyId, b3OffsetPos( far, (b3Vec3){ 0.0f, 0.0f, 1.0f } ), b3Quat_identity );

	for ( int i = 0; i < 600; ++i )
	{
		b3World_Step( worldId, 1.0f / 60.0f, 4 );
	}

	b3Pos p = b3Body_GetPosition( bodyId );
	ENSURE_SMALL( (float)( p.z - far.z ), 0.02f );
	ENSURE_SMALL( (float)( p.x - far.x ), 0.01f );

	b3DestroyWorld( worldId );
	return 0;
}
#endif

int FluidTest( void )
{
	RUN_SUBTEST( SubmergedVolumeTest );
	RUN_SUBTEST( BuoyancyTest );
	RUN_SUBTEST( FluidDragTest );
	RUN_SUBTEST( MultiShapeDragTest );
	RUN_SUBTEST( FluidWakeTest );
	RUN_SUBTEST( NoEventsTest );
	RUN_SUBTEST( LongRodTest );
#if defined( BOX3D_DOUBLE_PRECISION )
	RUN_SUBTEST( FarFluidTest );
#endif
	return 0;
}
