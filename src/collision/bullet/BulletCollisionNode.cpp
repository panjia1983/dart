/*
 * Copyright (c) 2011, Georgia Tech Research Corporation
 * All rights reserved.
 *
 * Author(s): Jeongseok Lee <jslee02@gmail.com>
 * Date: 05/11/2013
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

#include <assimp/scene.h>

#include "dynamics/BodyNode.h"
#include "dynamics/BoxShape.h"
#include "dynamics/EllipsoidShape.h"
#include "dynamics/CylinderShape.h"
#include "dynamics/PlaneShape.h"
#include "dynamics/MeshShape.h"

#include "collision/bullet/BulletCollisionNode.h"
#include "collision/bullet/BulletTypes.h"

namespace dart {
namespace collision {

BulletCollisionNode::BulletCollisionNode(dynamics::BodyNode* _bodyNode)
    : CollisionNode(_bodyNode)
{
    for(int i = 0; i < _bodyNode->getNumCollisionShapes(); i++)
    {
        dynamics::Shape* shape = _bodyNode->getCollisionShape(i);
        switch (shape->getShapeType())
        {
            case dynamics::Shape::BOX:
            {
                dynamics::BoxShape* box = dynamic_cast<dynamics::BoxShape*>(shape);

                btBoxShape* btBox = new btBoxShape(btVector3(box->getDim()[0]*0.5, box->getDim()[1]*0.5, box->getDim()[2]*0.5));
                btBox->setMargin(0.0);
                btCollisionObject* btCollObj = new btCollisionObject();
                btCollObj->setCollisionShape(btBox);
                btUserData* userData = new btUserData;
                userData->bodyNode = _bodyNode;
                userData->shape = shape;
                userData->btCollNode = this;
                btCollObj->setUserPointer(userData);
                mbtCollsionObjects.push_back(btCollObj);

                break;
            }
            case dynamics::Shape::ELLIPSOID:
            {
                dynamics::EllipsoidShape* ellipsoid = dynamic_cast<dynamics::EllipsoidShape*>(shape);

                if (ellipsoid->isSphere())
                {
                    btSphereShape* btCylinder = new btSphereShape(ellipsoid->getDim()[0] * 0.5);
                    btCylinder->setMargin(0.0);
                    btCollisionObject* btCollObj = new btCollisionObject();
                    btCollObj->setCollisionShape(btCylinder);
                    btUserData* userData = new btUserData;
                    userData->bodyNode = _bodyNode;
                    userData->shape = shape;
                    userData->btCollNode = this;
                    btCollObj->setUserPointer(userData);
                    mbtCollsionObjects.push_back(btCollObj);
                }
                else
                {
                }

                break;
            }
            case dynamics::Shape::CYLINDER:
            {
                dynamics::CylinderShape* cylinder = dynamic_cast<dynamics::CylinderShape*>(shape);

                btCylinderShapeZ* btCylinder = new btCylinderShapeZ(btVector3(cylinder->getRadius(), cylinder->getRadius(), cylinder->getHeight() * 0.5));
                btCylinder->setMargin(0.0);
                btCollisionObject* btCollObj = new btCollisionObject();
                btCollObj->setCollisionShape(btCylinder);
                btUserData* userData = new btUserData;
                userData->bodyNode = _bodyNode;
                userData->shape = shape;
                userData->btCollNode = this;
                btCollObj->setUserPointer(userData);
                mbtCollsionObjects.push_back(btCollObj);

                break;
            }
            case dynamics::Shape::PLANE:
            {
                dynamics::PlaneShape* plane = dynamic_cast<dynamics::PlaneShape*>(shape);

                double d = plane->getNormal().dot(plane->getPoint()) / plane->getNormal().squaredNorm();

                btStaticPlaneShape* btStaticPlane = new btStaticPlaneShape(convertVector3(plane->getNormal()), d);
                btStaticPlane->setMargin(0.0);
                btCollisionObject* btCollObj = new btCollisionObject();
                btCollObj->setCollisionShape(btStaticPlane);
                btUserData* userData = new btUserData;
                userData->bodyNode = _bodyNode;
                userData->shape = shape;
                userData->btCollNode = this;
                btCollObj->setUserPointer(userData);
                mbtCollsionObjects.push_back(btCollObj);

                break;
            }
            case dynamics::Shape::MESH:
            {
                dynamics::MeshShape *shapeMesh
                        = dynamic_cast<dynamics::MeshShape *>(shape);

                break;
            }
            default:
            {
                std::cout << "ERROR: Collision checking does not support "
                          << _bodyNode->getName()
                          << "'s Shape type\n";
                break;
            }
        }
    }
}

BulletCollisionNode::~BulletCollisionNode()
{
}

void BulletCollisionNode::updateBTCollisionObjects()
{
    for (int i = 0; i < mbtCollsionObjects.size(); ++i)
    {
        btUserData* userData = static_cast<btUserData*>(mbtCollsionObjects[i]->getUserPointer());
        dynamics::Shape* shape = userData->shape;
        btTransform T = convertTransform(mBodyNode->getWorldTransform() * shape->getLocalTransform());
        mbtCollsionObjects[i]->setWorldTransform(T);
    }
}

int BulletCollisionNode::getNumBTCollisionObjects() const
{
    return mbtCollsionObjects.size();
}

btCollisionObject*BulletCollisionNode::getBTCollisionObject(int _i)
{
    return mbtCollsionObjects[_i];
}

} // namespace collision
} // namespace dart
