// SPDX-FileCopyrightText: 2026 Erin Catto
// SPDX-License-Identifier: MIT

#include "body.h"
#include "core.h"
#include "joint.h"
#include "physics_world.h"
#include "recording.h"
#include "solver.h"
#include "solver_set.h"

// needed for dll export
#include "box3d/box3d.h"

// Pulley:
// length1 = norm(p1 - s1)
// length2 = norm(p2 - s2)
// C0 = (length1 + ratio * length2)_initial
// C = C0 - (length1 + ratio * length2)
// u1 = (p1 - s1) / norm(p1 - s1)
// u2 = (p2 - s2) / norm(p2 - s2)
// Cdot = -dot(u1, v1 + cross(w1, r1)) - ratio * dot(u2, v2 + cross(w2, r2))
// J = -[u1 cross(r1, u1) ratio * u2  ratio * cross(r2, u2)]
// K = J * invM * JT
//   = invMass1 + invI1 * cross(r1, u1)^2 + ratio^2 * (invMass2 + invI2 * cross(r2, u2)^2)

void b3PulleyJoint_SetLength( b3JointId jointId, float length )
{
	b3World* world = b3GetWorld( jointId.world0 );
	B3_REC( world, PulleyJointSetLength, jointId, length );
	b3JointSim* base = b3GetJointSimCheckType( jointId, b3_pulleyJoint );
	base->pulleyJoint.length = b3ClampFloat( length, B3_LINEAR_SLOP, B3_HUGE );
	base->pulleyJoint.impulse = 0.0f;
}

float b3PulleyJoint_GetLength( b3JointId jointId )
{
	b3JointSim* base = b3GetJointSimCheckType( jointId, b3_pulleyJoint );
	return base->pulleyJoint.length;
}

float b3PulleyJoint_GetRatio( b3JointId jointId )
{
	b3JointSim* base = b3GetJointSimCheckType( jointId, b3_pulleyJoint );
	return base->pulleyJoint.ratio;
}

b3Pos b3PulleyJoint_GetGroundAnchorA( b3JointId jointId )
{
	b3JointSim* base = b3GetJointSimCheckType( jointId, b3_pulleyJoint );
	return base->pulleyJoint.groundAnchorA;
}

b3Pos b3PulleyJoint_GetGroundAnchorB( b3JointId jointId )
{
	b3JointSim* base = b3GetJointSimCheckType( jointId, b3_pulleyJoint );
	return base->pulleyJoint.groundAnchorB;
}

float b3PulleyJoint_GetCurrentLength( b3JointId jointId )
{
	b3JointSim* base = b3GetJointSimCheckType( jointId, b3_pulleyJoint );

	b3World* world = b3GetUnlockedWorld( jointId.world0 );
	if ( world == NULL )
	{
		return 0.0f;
	}

	b3PulleyJoint* joint = &base->pulleyJoint;
	b3WorldTransform transformA = b3GetBodyTransform( world, base->bodyIdA );
	b3WorldTransform transformB = b3GetBodyTransform( world, base->bodyIdB );

	b3Pos pA = b3TransformWorldPoint( transformA, base->localFrameA.p );
	b3Pos pB = b3TransformWorldPoint( transformB, base->localFrameB.p );
	float lengthA = b3Length( b3SubPos( pA, joint->groundAnchorA ) );
	float lengthB = b3Length( b3SubPos( pB, joint->groundAnchorB ) );
	return lengthA + joint->ratio * lengthB;
}

b3Vec3 b3GetPulleyJointForce( b3World* world, b3JointSim* base )
{
	b3PulleyJoint* joint = &base->pulleyJoint;

	b3WorldTransform transformB = b3GetBodyTransform( world, base->bodyIdB );
	b3Pos pB = b3TransformWorldPoint( transformB, base->localFrameB.p );
	b3Vec3 uB = b3Normalize( b3SubPos( pB, joint->groundAnchorB ) );
	return b3MulSV( -joint->ratio * joint->impulse * world->inv_h, uB );
}

void b3PreparePulleyJoint( b3JointSim* base, b3StepContext* context )
{
	B3_ASSERT( base->type == b3_pulleyJoint );

	int idA = base->bodyIdA;
	int idB = base->bodyIdB;

	b3World* world = context->world;
	b3Body* bodyA = b3Array_Get( world->bodies, idA );
	b3Body* bodyB = b3Array_Get( world->bodies, idB );

	B3_ASSERT( bodyA->setIndex == b3_awakeSet || bodyB->setIndex == b3_awakeSet );

	b3SolverSet* setA = b3Array_Get( world->solverSets, bodyA->setIndex );
	b3SolverSet* setB = b3Array_Get( world->solverSets, bodyB->setIndex );

	int localIndexA = bodyA->localIndex;
	int localIndexB = bodyB->localIndex;

	b3BodySim* bodySimA = b3Array_Get( setA->bodySims, localIndexA );
	b3BodySim* bodySimB = b3Array_Get( setB->bodySims, localIndexB );

	float mA = bodySimA->invMass;
	b3Matrix3 iA = bodySimA->invInertiaWorld;
	float mB = bodySimB->invMass;
	b3Matrix3 iB = bodySimB->invInertiaWorld;

	base->invMassA = mA;
	base->invMassB = mB;
	base->invIA = iA;
	base->invIB = iB;

	b3PulleyJoint* joint = &base->pulleyJoint;

	joint->indexA = bodyA->setIndex == b3_awakeSet ? localIndexA : B3_NULL_INDEX;
	joint->indexB = bodyB->setIndex == b3_awakeSet ? localIndexB : B3_NULL_INDEX;

	joint->anchorA = b3RotateVector( bodySimA->transform.q, b3Sub( base->localFrameA.p, bodySimA->localCenter ) );
	joint->anchorB = b3RotateVector( bodySimB->transform.q, b3Sub( base->localFrameB.p, bodySimB->localCenter ) );
	joint->groundToCenterA = b3SubPos( bodySimA->center, joint->groundAnchorA );
	joint->groundToCenterB = b3SubPos( bodySimB->center, joint->groundAnchorB );

	b3Vec3 uA = b3Normalize( b3Add( joint->groundToCenterA, joint->anchorA ) );
	b3Vec3 uB = b3Normalize( b3Add( joint->groundToCenterB, joint->anchorB ) );

	// Compute effective mass.
	float ratio = joint->ratio;
	b3Vec3 crA = b3Cross( joint->anchorA, uA );
	b3Vec3 crB = b3Cross( joint->anchorB, uB );

	float kA = mA + b3Dot( crA, b3MulMV( iA, crA ) );
	float kB = mB + b3Dot( crB, b3MulMV( iB, crB ) );
	float k = kA + ratio * ratio * kB;
	joint->axialMass = k > 0.0f ? 1.0f / k : 0.0f;

	if ( context->enableWarmStarting == false )
	{
		joint->impulse = 0.0f;
	}
}

void b3WarmStartPulleyJoint( b3JointSim* base, b3StepContext* context )
{
	B3_ASSERT( base->type == b3_pulleyJoint );

	float mA = base->invMassA;
	float mB = base->invMassB;
	b3Matrix3 iA = base->invIA;
	b3Matrix3 iB = base->invIB;

	// Dummy state for static bodies
	b3BodyState dummyState = b3_identityBodyState;

	b3PulleyJoint* joint = &base->pulleyJoint;
	b3BodyState* stateA = joint->indexA == B3_NULL_INDEX ? &dummyState : context->states + joint->indexA;
	b3BodyState* stateB = joint->indexB == B3_NULL_INDEX ? &dummyState : context->states + joint->indexB;

	b3Vec3 rA = b3RotateVector( stateA->deltaRotation, joint->anchorA );
	b3Vec3 rB = b3RotateVector( stateB->deltaRotation, joint->anchorB );

	b3Vec3 uA = b3Normalize( b3Add( b3Add( joint->groundToCenterA, stateA->deltaPosition ), rA ) );
	b3Vec3 uB = b3Normalize( b3Add( b3Add( joint->groundToCenterB, stateB->deltaPosition ), rB ) );

	b3Vec3 PA = b3MulSV( -joint->impulse, uA );
	b3Vec3 PB = b3MulSV( -joint->ratio * joint->impulse, uB );

	if ( stateA->flags & b3_dynamicFlag )
	{
		stateA->linearVelocity = b3MulAdd( stateA->linearVelocity, mA, PA );
		stateA->angularVelocity = b3Add( stateA->angularVelocity, b3MulMV( iA, b3Cross( rA, PA ) ) );
	}

	if ( stateB->flags & b3_dynamicFlag )
	{
		stateB->linearVelocity = b3MulAdd( stateB->linearVelocity, mB, PB );
		stateB->angularVelocity = b3Add( stateB->angularVelocity, b3MulMV( iB, b3Cross( rB, PB ) ) );
	}
}

void b3SolvePulleyJoint( b3JointSim* base, b3StepContext* context, bool useBias )
{
	B3_ASSERT( base->type == b3_pulleyJoint );

	float mA = base->invMassA;
	float mB = base->invMassB;
	b3Matrix3 iA = base->invIA;
	b3Matrix3 iB = base->invIB;

	// Dummy state for static bodies
	b3BodyState dummyState = b3_identityBodyState;

	b3PulleyJoint* joint = &base->pulleyJoint;
	b3BodyState* stateA = joint->indexA == B3_NULL_INDEX ? &dummyState : context->states + joint->indexA;
	b3BodyState* stateB = joint->indexB == B3_NULL_INDEX ? &dummyState : context->states + joint->indexB;

	b3Vec3 vA = stateA->linearVelocity;
	b3Vec3 wA = stateA->angularVelocity;
	b3Vec3 vB = stateB->linearVelocity;
	b3Vec3 wB = stateB->angularVelocity;

	// anchors
	b3Vec3 rA = b3RotateVector( stateA->deltaRotation, joint->anchorA );
	b3Vec3 rB = b3RotateVector( stateB->deltaRotation, joint->anchorB );

	// rope segments
	b3Vec3 dA = b3Add( b3Add( joint->groundToCenterA, stateA->deltaPosition ), rA );
	b3Vec3 dB = b3Add( b3Add( joint->groundToCenterB, stateB->deltaPosition ), rB );

	float lengthA, lengthB;
	b3Vec3 uA = b3GetLengthAndNormalize( &lengthA, dA );
	b3Vec3 uB = b3GetLengthAndNormalize( &lengthB, dB );

	float ratio = joint->ratio;
	float C = joint->length - lengthA - ratio * lengthB;

	b3Vec3 vpA = b3Add( vA, b3Cross( wA, rA ) );
	b3Vec3 vpB = b3Add( vB, b3Cross( wB, rB ) );
	float Cdot = -b3Dot( uA, vpA ) - ratio * b3Dot( uB, vpB );

	float bias = 0.0f;
	float massScale = 1.0f;
	float impulseScale = 0.0f;
	bool slack = joint->isRigid == false && C > ( joint->holdTaut ? B3_LINEAR_SLOP : 0.0f );
	if ( slack )
	{
		bias = C * context->inv_h;
	}
	else if ( useBias )
	{
		bias = base->constraintSoftness.biasRate * C;
		massScale = base->constraintSoftness.massScale;
		impulseScale = base->constraintSoftness.impulseScale;
	}

	float impulse = -massScale * joint->axialMass * ( Cdot + bias ) - impulseScale * joint->impulse;
	float oldImpulse = joint->impulse;
	joint->impulse = slack ? b3MaxFloat( 0.0f, oldImpulse + impulse ) : oldImpulse + impulse;
	impulse = joint->impulse - oldImpulse;

	b3Vec3 PA = b3MulSV( -impulse, uA );
	b3Vec3 PB = b3MulSV( -ratio * impulse, uB );
	vA = b3MulAdd( vA, mA, PA );
	wA = b3Add( wA, b3MulMV( iA, b3Cross( rA, PA ) ) );
	vB = b3MulAdd( vB, mB, PB );
	wB = b3Add( wB, b3MulMV( iB, b3Cross( rB, PB ) ) );

	if ( stateA->flags & b3_dynamicFlag )
	{
		stateA->linearVelocity = vA;
		stateA->angularVelocity = wA;
	}

	if ( stateB->flags & b3_dynamicFlag )
	{
		stateB->linearVelocity = vB;
		stateB->angularVelocity = wB;
	}
}

void b3DrawPulleyJoint( b3DebugDraw* draw, b3JointSim* base, b3WorldTransform transformA, b3WorldTransform transformB )
{
	B3_ASSERT( base->type == b3_pulleyJoint );

	b3PulleyJoint* joint = &base->pulleyJoint;

	b3Pos pA = b3TransformWorldPoint( transformA, base->localFrameA.p );
	b3Pos pB = b3TransformWorldPoint( transformB, base->localFrameB.p );

	draw->DrawSegmentFcn( joint->groundAnchorA, pA, b3_colorWhite, draw->context );
	draw->DrawSegmentFcn( joint->groundAnchorA, joint->groundAnchorB, b3_colorGray, draw->context );
	draw->DrawSegmentFcn( joint->groundAnchorB, pB, b3_colorWhite, draw->context );
	draw->DrawPointFcn( joint->groundAnchorA, 6.0f, b3_colorLightGreen, draw->context );
	draw->DrawPointFcn( joint->groundAnchorB, 6.0f, b3_colorLightGreen, draw->context );
	draw->DrawPointFcn( pA, 4.0f, b3_colorWhite, draw->context );
	draw->DrawPointFcn( pB, 4.0f, b3_colorWhite, draw->context );
}
