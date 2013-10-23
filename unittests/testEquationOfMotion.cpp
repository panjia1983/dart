/*
 * Copyright (c) 2011, Georgia Tech Research Corporation
 * All rights reserved.
 *
 * Author(s): Jeongseok Lee <jslee02@gmail.com>
 * Date: 05/23/2013
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
#include "dynamics/FreeJoint.h"
#include "dynamics/PrismaticJoint.h"
#include "dynamics/RevoluteJoint.h"
#include "dynamics/Skeleton.h"
#include "dynamics/TranslationalJoint.h"
#include "dynamics/UniversalJoint.h"
#include "dynamics/WeldJoint.h"
#include "dynamics/EulerJoint.h"
#include "dynamics/ScrewJoint.h"
#include "simulation/World.h"
#include "utils/Paths.h"
#include "utils/SkelParser.h"

using namespace dart;
using namespace math;
using namespace dynamics;

#define EOM_TOL 0.01

/******************************************************************************/
TEST(EOM, SinglePendulum)
{
    simulation::World* myWorld = utils::SkelParser::readSkelFile(
            DART_DATA_PATH"/skel/test/single_pendulum.skel");
    EXPECT_TRUE(myWorld != NULL);

    myWorld->setGravity(Eigen::Vector3d(0.0, 0.0, -9.81));

    dynamics::Skeleton* pendulum = myWorld->getSkeleton("single_pendulum");
    EXPECT_TRUE(pendulum != NULL);

    double simTime = 2.0;
    double timeStep = myWorld->getTimeStep();
    int nSteps = simTime / timeStep;

    for (int i = 0; i < nSteps; i++)
    {
        myWorld->step();

        // mass matrix
        Eigen::MatrixXd M = pendulum->getMassMatrix();
        Eigen::MatrixXd M_FS = pendulum->getMassMatrixFS();

        if (!equals(M, M_FS))
        {
            std::cout << "M   : " << M << std::endl;
            std::cout << "M_FS: " << M_FS << std::endl;
        }

        EXPECT_TRUE(equals(M, M_FS));

        // inverse mass matrix
//        Eigen::MatrixXd MInv = pendulum->getInvMassMatrix();
//        Eigen::MatrixXd MInv_FS = pendulum->getInvMassMatrixFS();

//        if (!equals(MInv, MInv_FS))
//        {
//            std::cout << "MInv   : \n" << MInv << std::endl;
//            std::cout << "MInv_FS: \n" << MInv_FS << std::endl;
//        }

//        EXPECT_TRUE(equals(MInv, MInv_FS));
    }
}

/******************************************************************************/
TEST(EOM, DoulbePendulum)
{
    double tol = 1e-4;

    simulation::World* myWorld = utils::SkelParser::readSkelFile(
            DART_DATA_PATH"/skel/test/double_pendulum.skel");
    EXPECT_TRUE(myWorld != NULL);

    myWorld->setGravity(Eigen::Vector3d(0.0, 0.0, -9.81));

    dynamics::Skeleton* pendulum = myWorld->getSkeleton("double_pendulum");
    EXPECT_TRUE(pendulum != NULL);

    double simTime = 0.001;
    double timeStep = myWorld->getTimeStep();
    int nSteps = simTime / timeStep;

    for (int i = 0; i < nSteps; i++)
    {
        myWorld->step();

        // mass matrix
        Eigen::MatrixXd M = pendulum->getMassMatrix();
        Eigen::MatrixXd M_FS = pendulum->getMassMatrixFS();

        if (!equals(M, M_FS))
        {
            std::cout << "M   : \n" << M << std::endl;
            std::cout << "M_FS: \n" << M_FS << std::endl;
        }

        EXPECT_TRUE(equals(M, M_FS));

        // inverse mass matrix
//        Eigen::MatrixXd MInv = pendulum->getInvMassMatrix();
//        Eigen::MatrixXd MInv_FS = pendulum->getInvMassMatrixFS();

//        if (!equals(MInv, MInv_FS))
//        {
//            std::cout << "MInv   : \n" << MInv << std::endl;
//            std::cout << "MInv_FS: \n" << MInv_FS << std::endl;
//        }

//        EXPECT_TRUE(equals(MInv, MInv_FS));
    }
}

/******************************************************************************/
TEST(EOM, FullBody1)
{
    double tol = 1e-4;

    simulation::World* myWorld = utils::SkelParser::readSkelFile(
            DART_DATA_PATH"/skel/fullbody1.skel");
    EXPECT_TRUE(myWorld != NULL);

    myWorld->setGravity(Eigen::Vector3d(0.0, 0.0, -9.81));

    dynamics::Skeleton* pendulum = myWorld->getSkeleton(0);
    EXPECT_TRUE(pendulum != NULL);

    double simTime = 2.000;
    double timeStep = myWorld->getTimeStep();
    int nSteps = simTime / timeStep;

    for (int i = 0; i < nSteps; i++)
    {
        myWorld->step();

        Eigen::MatrixXd M = pendulum->getMassMatrix();
        Eigen::MatrixXd M_FS = pendulum->getMassMatrixFS();

        if (!equals(M, M_FS))
        {
            std::cout << "M   : \n" << M << std::endl;
            std::cout << "M_FS: \n" << M_FS << std::endl;
        }

        EXPECT_TRUE(equals(M, M_FS));

//        Eigen::MatrixXd MInv = pendulum->getInvMassMatrix();
//        Eigen::MatrixXd MInv_FS = pendulum->getInvMassMatrixFS();

//        if (!equals(MInv, MInv_FS))
//        {
//            std::cout << "MInv   : \n" << MInv << std::endl;
//            std::cout << "MInv_FS: \n" << MInv_FS << std::endl;
//        }

//        EXPECT_TRUE(equals(MInv, MInv_FS));
    }
}

/******************************************************************************/
int main(int argc, char* argv[])
{
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}


