/*
 * Copyright (c) 2011, Georgia Tech Research Corporation
 * All rights reserved.
 *
 * Author(s): Jeongseok Lee <jslee02@gmail.com>
 * Date: 10/11/2013
 *
 * Geoorgia Tech Graphics Lab and Humanoid Robotics Lab
 *
 * Directed by Prof. C. Karen Liu and Prof. Mike Stilman
 * <karenliu@cc.gatech.edu> <mstilman@cc.gatech.edu>
 *
 * This file is provided under the following "BSD-style" License:
 *   Redistribution and use in source and binary forms, with or
 *   without modification, are permitted provided that the following
 *   conditions are met:
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
 *   CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
 *   INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *   MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 *   DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
 *   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 *   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
 *   USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *   AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *   LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *   ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *   POSSIBILITY OF SUCH DAMAGE.
 */

#include <iostream>
#include <gtest/gtest.h>
#include "TestHelpers.h"

#include "math/Geometry.h"
#include "math/Helpers.h"
#include "dynamics/BallJoint.h"
#include "dynamics/RevoluteJoint.h"
#include "dynamics/PrismaticJoint.h"
#include "dynamics/FreeJoint.h"
#include "dynamics/WeldJoint.h"
#include "dynamics/TranslationalJoint.h"
#include "dynamics/Skeleton.h"
#include "simulation/World.h"
#include "utils/Paths.h"
#include "utils/SkelParser.h"

using namespace dart;
using namespace math;
using namespace dynamics;
using namespace simulation;
using namespace constraint;
using namespace collision;
using namespace utils;

#define IMPULSE_TOL 1e-6

/******************************************************************************/
TEST(IMPULSE, BASIC)
{
    int numTest = 100;
    double desiredVel0 = -10;

    World* world = SkelParser::readSkelFile(DART_DATA_PATH"skel/test/impulse_test.skel");
    EXPECT_TRUE(world != NULL);

    double timeStep = world->getTimeStep();

    Eigen::Vector3d gravity = world->getGravity();

    Skeleton* sphere1Skeleton = world->getSkeleton(1);
    EXPECT_TRUE(sphere1Skeleton != NULL);

    BodyNode* sphere1BodyNode = sphere1Skeleton->getBodyNode(0);
    EXPECT_TRUE(sphere1BodyNode != NULL);

    Eigen::VectorXd state = sphere1Skeleton->getState();
    EXPECT_TRUE(state.size() == 12);
    state[11] = desiredVel0;
    sphere1Skeleton->setState(state);

    double mass = sphere1BodyNode->getMass();

    double posZ0 = sphere1BodyNode->getWorldTransform().translation()[2];
    EXPECT_EQ(posZ0, 0.5);

    double velZ0 = sphere1BodyNode->getWorldVelocity()[5];
    EXPECT_EQ(velZ0, desiredVel0);

    double accZ0 = sphere1BodyNode->getWorldAcceleration()[5];
    EXPECT_EQ(accZ0, 0.0);

    //--------------------------------------------------------------------------
    world->step();

    ConstraintDynamics* constraintDynamics = world->getConstraintHandler();
    EXPECT_TRUE(constraintDynamics != NULL);

    CollisionDetector* collisionDetector = constraintDynamics->getCollisionDetector();
    EXPECT_TRUE(collisionDetector != NULL);

    int numContact = collisionDetector->getNumContacts();
    EXPECT_EQ(numContact, 1);

    Contact contact = collisionDetector->getContact(0);

    double penetrationDepth = contact.penetrationDepth;
    EXPECT_EQ(penetrationDepth, 0.0);

    Eigen::Vector3d contactPoint = contact.point;
    EXPECT_EQ(contactPoint, Eigen::Vector3d::Zero());

    Eigen::Vector3d contactForce = contact.force;
    // We need 1e-1 tolerance because we added small positive value to A matrix
    // of LCP.
    // see: ConstraintDynamics::fillMatricesODE()
    EXPECT_NEAR(contactForce[2], mass * gravity[2] + mass * desiredVel0 / timeStep, 1e-6);

    double posZ = sphere1BodyNode->getWorldTransform().translation()[2];
    EXPECT_NEAR(posZ, 0.5, 1e-6);

    double velZ = sphere1BodyNode->getWorldVelocity()[5];
    EXPECT_NEAR(velZ, 0.0, 1e-6);

    double accZ = sphere1BodyNode->getWorldAcceleration()[5];
    EXPECT_NEAR(accZ, (mass * gravity[2] - contactForce[2]) / mass, 1e-6);

    //--------------------------------------------------------------------------
    world->step();

    numContact = collisionDetector->getNumContacts();
    EXPECT_EQ(numContact, 1);

    contact = collisionDetector->getContact(0);

    penetrationDepth = contact.penetrationDepth;
    EXPECT_EQ(penetrationDepth, 0.0);

    contactPoint = contact.point;
    EXPECT_EQ(contactPoint, Eigen::Vector3d::Zero());

    contactForce = contact.force;
    // We need 1e-1 tolerance because we added small positive value to A matrix
    // of LCP.
    // see: ConstraintDynamics::fillMatricesODE()
    EXPECT_NEAR(contactForce[2], mass * gravity[2], 1e-6);

    posZ = sphere1BodyNode->getWorldTransform().translation()[2];
    EXPECT_NEAR(posZ, 0.5, 1e-6);

    velZ = sphere1BodyNode->getWorldVelocity()[5];
    EXPECT_NEAR(velZ, 0.0, 1e-6);

    accZ = sphere1BodyNode->getWorldAcceleration()[5];
    EXPECT_NEAR(accZ, (mass * gravity[2] - contactForce[2]) / mass, 1e-6);

    //--------------------------------------------------------------------------
    world->step();

    numContact = collisionDetector->getNumContacts();
    EXPECT_EQ(numContact, 1);

    contact = collisionDetector->getContact(0);

    penetrationDepth = contact.penetrationDepth;
    EXPECT_EQ(penetrationDepth, 0.0);

    contactPoint = contact.point;
    EXPECT_EQ(contactPoint, Eigen::Vector3d::Zero());

    contactForce = contact.force;
    // We need 1e-1 tolerance because we added small positive value to A matrix
    // of LCP.
    // see: ConstraintDynamics::fillMatricesODE()
    EXPECT_NEAR(contactForce[2], mass * gravity[2], 1e-6);

    posZ = sphere1BodyNode->getWorldTransform().translation()[2];
    EXPECT_NEAR(posZ, 0.5, 1e-6);

    velZ = sphere1BodyNode->getWorldVelocity()[5];
    EXPECT_NEAR(velZ, 0.0, 1e-6);

    accZ = sphere1BodyNode->getWorldAcceleration()[5];
    EXPECT_NEAR(accZ, (mass * gravity[2] - contactForce[2]) / mass, 1e-6);
}

/******************************************************************************/
int main(int argc, char* argv[])
{
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}


