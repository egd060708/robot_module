#pragma once

#include <Eigen/Dense>
#include "math.h"
#include <iostream>
#include "SecondButterworthLPF.h"

using namespace Eigen;

namespace Quadruped
{
    // 单腿数据结构
    typedef struct _LegS
    {
        Vector3d Position = Vector3d::Zero();
        Vector3d Velocity = Vector3d::Zero();// 雅可比矩阵映射速度
        Vector3d Acc = Vector3d::Zero();
        Vector3d VelocityW = Vector3d::Zero();// 旋转速度
        Vector3d VelocityG = Vector3d::Zero();// 合速度
        Vector3d Force = Vector3d::Zero();// 执行的力
        Vector3d extForce = Vector3d::Zero();// 与外界的接触力

    } LegS;
    // 关节数据结构
    typedef struct _JointS
    {
        // 三轴单腿结构
        Vector3d Angle = Vector3d::Zero();
        Vector3d Velocity = Vector3d::Zero();
        Vector3d Acc = Vector3d::Zero();
        Vector3d Torque = Vector3d::Zero();
        // 末端执行器结构
        double Foot_Angle = 0;
        double Foot_Velocity = 0;
        double Foot_Acc = 0;
        double Foot_Torque = 0;
    } JointS;

    // 单腿类型
    class Leg
    {
    public:
        LegS targetLeg;
        LegS currentLeg;
        JointS targetJoint;
        JointS currentJoint;
        Matrix3d jacobi = Matrix3d::Zero();
        Matrix3d jacobi_dot = Matrix3d::Zero();
        double L1, L2, L2b, L3, L3b, L4, L4b;// 前三条连杆以及最后两个足端执行器参数
        double Reff; // 末端执行器触地点速度计算等效半径
        Vector3d Pcj[4];// 所有连杆的质心相对其关节位置
        Vector3d Pcl[4];// 所有连杆的质心相对单腿系的位置
        Vector<double, 6> Ic[4];// 所有连杆质心的惯性张量矩阵简化形式(基于连杆原始坐标系)
        Vector<double, 6> Icleg[4];// 所有连杆质心的惯性张量矩阵简化形式(基于单腿坐标系，其中单腿坐标系与机身坐标系平行)
        double Mc[4];// 所有连杆质量
        LPF_SecondOrder_Classdef enVelF[4] = { LPF_SecondOrder_Classdef(50,500),LPF_SecondOrder_Classdef(50,500),LPF_SecondOrder_Classdef(50,500),LPF_SecondOrder_Classdef(50,500) };

    public:
        Leg(double _l[7], Vector3d _pcj[4], double _mc[4], Vector<double, 6> _ic[4], int idx) {
            if (idx == 0 || idx == 2)
            {
                L1 = _l[0];
                L2b = _l[2];
                L3b = _l[4];
            }
            else
            {
                L1 = -_l[0];
                L2b = -_l[2];
                L3b = -_l[4];
            }
            L2 = _l[1];
            L3 = _l[3];
            L4 = _l[5];
            L4b = _l[6];
            for (int i = 0; i < 4; i++)
            {
                this->Mc[i] = _mc[i];
                this->Ic[i] = _ic[i];
                this->Icleg[i] = _ic[i];
                if (i == 0)
                {
                    switch (idx)
                    {
                    case 0:
                        Pcj[i](0) = _pcj[i](0);
                        Pcj[i](1) = _pcj[i](1);
                        Pcj[i](2) = _pcj[i](2);
                        break;
                    case 1:
                        Pcj[i](0) = _pcj[i](0);
                        Pcj[i](1) = -_pcj[i](1);
                        Pcj[i](2) = _pcj[i](2);
                        Ic[i](3) = -Ic[i](3);
                        Ic[i](5) = -Ic[i](5);
                        break;
                    case 2:
                        Pcj[i](0) = -_pcj[i](0);
                        Pcj[i](1) = _pcj[i](1);
                        Pcj[i](2) = _pcj[i](2);
                        Ic[i](3) = -Ic[i](3);
                        Ic[i](4) = -Ic[i](4);
                        break;
                    case 3:
                        Pcj[i](0) = -_pcj[i](0);
                        Pcj[i](1) = -_pcj[i](1);
                        Pcj[i](2) = _pcj[i](2);
                        Ic[i](4) = -Ic[i](4);
                        Ic[i](5) = -Ic[i](5);
                        break;
                    default:
                        break;
                    }
                }
                else
                {
                    switch (idx)
                    {
                    case 0:
                        Pcj[i](0) = _pcj[i](0);
                        Pcj[i](1) = _pcj[i](1);
                        Pcj[i](2) = _pcj[i](2);
                        break;
                    case 1:
                        Pcj[i](0) = _pcj[i](0);
                        Pcj[i](1) = -_pcj[i](1);
                        Pcj[i](2) = _pcj[i](2);
                        Ic[i](3) = -Ic[i](3);
                        Ic[i](5) = -Ic[i](5);
                        break;
                    case 2:
                        Pcj[i](0) = _pcj[i](0);
                        Pcj[i](1) = _pcj[i](1);
                        Pcj[i](2) = _pcj[i](2);
                        break;
                    case 3:
                        Pcj[i](0) = _pcj[i](0);
                        Pcj[i](1) = -_pcj[i](1);
                        Pcj[i](2) = _pcj[i](2);
                        Ic[i](3) = -Ic[i](3);
                        Ic[i](5) = -Ic[i](5);
                        break;
                    default:
                        break;
                    }
                }
            }
            Reff = L4 + L4b;
        }
        void updateJointAng(Vector4d _jAngle);          // 更新各关节观测角度
        void updateJointVel(Vector4d _jVel);            // 更新关节观测角速度
        void updateJointAcc(Vector4d _jAcc);            // 更新关节观测角加速度
        void updateJointTau(Vector4d _jTorque);         // 更新关节观测力矩
        void updateJointAng(Vector3d _jAngle);          // 更新各关节观测角度
        void updateJointVel(Vector3d _jVel);            // 更新关节观测角速度
        void updateJointAcc(Vector3d _jAcc);            // 更新关节观测角角加速度
        void updateJointTau(Vector3d _jTorque);         // 更新关节观测力矩
        void legJacobi_Cal();                           // 根据当前电机角速度计算雅可比矩阵及其微分矩阵
        void legFK_Cal();                               //不涉及接触点牵连速度的正运动学映射
        void legFK_Cal(const Vector3d& _imu, const Vector3d& _gyro);// 包括由当前角度映射末端位姿，由当前角速度映射末端速度，由当前力矩映射当前末端力
        void setTargetLegPositon(Vector3d _lPosition);  // 更新腿部末端位置目标值
        void setTargetLegVelocity(Vector3d _lVelocity); // 更新腿部末端速度目标值
        void setTargetLegForce(Vector3d _lForce);       // 更新腿部末端力目标值
        void setTargetEndTau(double _eTau);             // 更新末端执行器目标力矩
        void legIK_Cal();                               // 包括由末端目标位姿映射关节目标角度，由末端目标速度映射到关节目标速度，由末端目标力映射到当前关节力矩
        void legIP_Cal();                               // 腿部连杆惯量动态位置计算(从初始位型到当前位型)

        const JointS& getJointCurrent();                // 获取关节参数观测值
        const LegS& getLegCurrent();                    // 获取单腿模型参数观测值

        Eigen::Matrix3d aI2mI(const Vector<double, 6>& oriIa)
        {
            Matrix3d mI;
            mI(0, 0) = oriIa(0);
            mI(1, 1) = oriIa(1);
            mI(2, 2) = oriIa(2);
            mI(0, 1) = oriIa(3);
            mI(1, 0) = oriIa(3);
            mI(0, 2) = oriIa(4);
            mI(2, 0) = oriIa(4);
            mI(1, 2) = oriIa(5);
            mI(2, 1) = oriIa(5);
            return mI;
        }

        Eigen::Vector<double, 6> mI2aI(const Matrix3d& oriIm)
        {
            Vector<double, 6> aI;
            aI(0) = oriIm(0, 0);
            aI(1) = oriIm(1, 1);
            aI(2) = oriIm(2, 2);
            aI(3) = oriIm(0, 1);
            aI(4) = oriIm(0, 2);
            aI(5) = oriIm(1, 2);
            return aI;
        }
    };

    void Leg::updateJointAng(Vector4d _jAngle)
    {
        currentJoint.Angle = _jAngle.block(0, 0, 3, 1);
        currentJoint.Foot_Angle = _jAngle(3);
    }

    void Leg::updateJointVel(Vector4d _jVel)
    {
        currentJoint.Velocity = _jVel.block(0, 0, 3, 1);
        currentJoint.Foot_Velocity = _jVel(3);
    }

    void Leg::updateJointAcc(Vector4d _jAcc)
    {
        currentJoint.Acc = _jAcc.block(0, 0, 3, 1);
        currentJoint.Foot_Acc = _jAcc(3);
    }

    void Leg::updateJointTau(Vector4d _jTorque)
    {
        currentJoint.Torque = _jTorque.block(0,0,3,1);
        currentJoint.Foot_Torque = _jTorque(3);
    }

    void Leg::updateJointAng(Vector3d _jAngle)
    {
        currentJoint.Angle = _jAngle;
    }

    void Leg::updateJointVel(Vector3d _jVel)
    {
        currentJoint.Velocity = _jVel;
    }

    void Leg::updateJointAcc(Vector3d _jAcc)
    {
        currentJoint.Acc = _jAcc;
    }

    void Leg::updateJointTau(Vector3d _jTorque)
    {
        currentJoint.Torque = _jTorque;
    }

    void Leg::legJacobi_Cal()
    {
        // 雅可比矩阵计算
        Matrix3d J = Matrix3d::Zero();
        Matrix3d J_dot = Matrix3d::Zero();
        double s1 = sin(currentJoint.Angle(0));
        double c1 = cos(currentJoint.Angle(0));
        double s2 = sin(currentJoint.Angle(1));
        double c2 = cos(currentJoint.Angle(1));
        double s23 = sin(currentJoint.Angle(1) + currentJoint.Angle(2));
        double c23 = cos(currentJoint.Angle(1) + currentJoint.Angle(2));
        J(0, 0) = 0;
        J(0, 1) = -L2 * c2 - L3 * c23;
        J(0, 2) = -L3 * c23;
        J(1, 0) = -L1 * s1 + L2 * c1 * c2 + L3 * c1 * c23 - L2b * s1 - L3b * s1 + L4 * c1;
        J(1, 1) = -L2 * s1 * s2 - L3 * s1 * s23;
        J(1, 2) = -L3 * s1 * s23;
        J(2, 0) = L1 * c1 + L2 * s1 * c2 + L3 * s1 * c23 + L2b * c1 + L3b * c1 + L4 * s1;
        J(2, 1) = L2 * c1 * s2 + L3 * c1 * s23;
        J(2, 2) = L3 * c1 * s23;
        this->jacobi = J;
        double s1_dot = c1 * currentJoint.Velocity(0);
        double c1_dot = -s1 * currentJoint.Velocity(0);
        double s2_dot = c2 * currentJoint.Velocity(1);
        double c2_dot = -s2 * currentJoint.Velocity(1);
        double s23_dot = c23 * (currentJoint.Velocity(1) + currentJoint.Velocity(2));
        double c23_dot = -s23 * (currentJoint.Velocity(1) + currentJoint.Velocity(2));
        J_dot(0, 0) = 0;
        J_dot(0, 1) = -L2 * c2_dot - L3 * c23_dot;
        J_dot(0, 2) = -L3 * c23_dot;
        J_dot(1, 0) = -L1 * s1_dot + L2 * c1_dot * c2_dot + L3 * c1_dot * c23_dot - L2b * s1_dot - L3b * s1_dot + L4 * c1_dot;
        J_dot(1, 1) = -L2 * s1_dot * s2_dot - L3 * s1_dot * s23_dot;
        J_dot(1, 2) = -L3 * s1_dot * s23_dot;
        J_dot(2, 0) = L1 * c1_dot + L2 * s1_dot * c2_dot + L3 * s1_dot * c23_dot + L2b * c1_dot + L3b * c1_dot + L4 * s1_dot;
        J_dot(2, 1) = L2 * c1_dot * s2_dot + L3 * c1_dot * s23_dot;
        J_dot(2, 2) = L3 * c1_dot * s23_dot;
        this->jacobi_dot = J_dot;
    }

    void Leg::legFK_Cal()
    {
        // 正运动学更新当前末端位置
        currentLeg.Position(0) = -L2 * sin(currentJoint.Angle(1)) - L3 * sin(currentJoint.Angle(1) + currentJoint.Angle(2));
        currentLeg.Position(1) = L1 * cos(currentJoint.Angle(0)) + L2 * sin(currentJoint.Angle(0)) * cos(currentJoint.Angle(1)) + L2b * cos(currentJoint.Angle(0)) + L3 * sin(currentJoint.Angle(0)) * cos(currentJoint.Angle(1) + currentJoint.Angle(2)) + L3b * cos(currentJoint.Angle(0)) + L4 * sin(currentJoint.Angle(0));
        currentLeg.Position(2) = L1 * sin(currentJoint.Angle(0)) - L2 * cos(currentJoint.Angle(0)) * cos(currentJoint.Angle(1)) + L2b * sin(currentJoint.Angle(0)) - L3 * cos(currentJoint.Angle(0)) * cos(currentJoint.Angle(1) + currentJoint.Angle(2)) + L3b * sin(currentJoint.Angle(0)) - L4 * cos(currentJoint.Angle(0)) - L4b;
        // 雅可比矩阵完成从关节速度到末端速度的映射
        currentLeg.Velocity = jacobi * currentJoint.Velocity;
        // 雅可比矩阵完成从当前关节力矩到当前末端虚拟力的映射
        currentLeg.Force = jacobi.transpose().inverse() * currentJoint.Torque;
        // 得到仅由关节贡献的腿部末端加速度
        currentLeg.Acc = jacobi_dot * currentJoint.Velocity + jacobi * currentJoint.Acc;
    }

    void Leg::legFK_Cal(const Vector3d& _imu, const Vector3d& _gyro)
    {
        // 正运动学更新当前末端位置
        currentLeg.Position(0) = -L2 * sin(currentJoint.Angle(1)) - L3 * sin(currentJoint.Angle(1) + currentJoint.Angle(2));
        currentLeg.Position(1) = L1 * cos(currentJoint.Angle(0)) + L2 * sin(currentJoint.Angle(0)) * cos(currentJoint.Angle(1)) + L2b * cos(currentJoint.Angle(0)) + L3 * sin(currentJoint.Angle(0)) * cos(currentJoint.Angle(1) + currentJoint.Angle(2)) + L3b * cos(currentJoint.Angle(0)) + L4 * sin(currentJoint.Angle(0));
        currentLeg.Position(2) = L1 * sin(currentJoint.Angle(0)) - L2 * cos(currentJoint.Angle(0)) * cos(currentJoint.Angle(1)) + L2b * sin(currentJoint.Angle(0)) - L3 * cos(currentJoint.Angle(0)) * cos(currentJoint.Angle(1) + currentJoint.Angle(2)) + L3b * sin(currentJoint.Angle(0)) - L4 * cos(currentJoint.Angle(0)) - L4b;
        // 更新末端执行器等效半径
        Reff = L4 + L4b - L4b * sin(currentJoint.Angle(0) + _imu(0));
        // 雅可比矩阵完成从关节速度到末端速度的映射
        currentLeg.Velocity = jacobi * currentJoint.Velocity;
        // 接触点转动牵连速度
        currentLeg.VelocityW << Reff * (_gyro(1) * cos(currentJoint.Angle(0)) + currentJoint.Velocity(1) + currentJoint.Velocity(2) + currentJoint.Foot_Velocity), L4b* (_gyro(0) + currentJoint.Velocity(0)), 0;
        //currentLeg.VelocityW << Reff * (currentJoint.Foot_Velocity), 0, 0;
        //std::cout << "currentLeg.Vel: \n" << currentLeg.VelocityW << std::endl;
        // 得到总的接触速度
        currentLeg.VelocityG = currentLeg.Velocity + currentLeg.VelocityW;
        // 雅可比矩阵完成从当前关节力矩到当前末端虚拟力的映射
        currentLeg.Force = jacobi.transpose().inverse() * currentJoint.Torque;
        // 得到仅由关节贡献的腿部末端加速度
        currentLeg.Acc = jacobi_dot * currentJoint.Velocity + jacobi * currentJoint.Acc;
    }

    void Leg::setTargetLegPositon(Vector3d _lPosition)
    {
        targetLeg.Position = _lPosition;
    }

    void Leg::setTargetLegVelocity(Vector3d _lVelocity)
    {
        targetLeg.Velocity = _lVelocity;
    }

    void Leg::setTargetLegForce(Vector3d _lForce)
    {
        targetLeg.Force = _lForce;
        //targetJoint.Torque = jacobi.inverse() * targetLeg.Force;
        targetJoint.Torque = jacobi.transpose() * targetLeg.Force;
    }

    void Leg::setTargetEndTau(double _eTau)
    {
        targetJoint.Foot_Torque = _eTau;
    }

    void Leg::legIK_Cal()
    {
        // 解析几何法逆运动学
        double x = targetLeg.Position(0);
        double y = targetLeg.Position(1) - L2b * cos(currentJoint.Angle(0)) - L3b * cos(currentJoint.Angle(0)) - L4 * sin(currentJoint.Angle(0));
        double z = targetLeg.Position(2) - L2b * sin(currentJoint.Angle(0)) - L3b * sin(currentJoint.Angle(0)) + L4 * cos(currentJoint.Angle(0)) + L4b;

        double L = sqrt(y * y + z * z - L1 * L1);// 首先是二三级连杆投影到yz平面的映射长度
        double theta1 = atan2((L1 * z + L * y), -L * z + L1 * y);
        L = sqrt(x * x + y * y + z * z - L1 * L1);// 随后是投影到三维空间的等效长度
        double theta3 = -3.1415926 + acos((L2 * L2 + L3 * L3 - L * L) / 2 / L2 / L3);

        double a1 = y * sin(theta1) - z * cos(theta1);
        double m1 = -L3 * sin(theta3);
        double m2 = -L3 * cos(theta3) - L2;
        double theta2 = atan2(a1 * m1 + x * m2, x * m1 - a1 * m2);

        targetJoint.Angle(0) = theta1;
        targetJoint.Angle(1) = theta2;
        targetJoint.Angle(2) = theta3;
    }

    void Leg::legIP_Cal()
    {
        AngleAxisd rhip(currentJoint.Angle(0), Vector3d::UnitX());
        AngleAxisd rthigh(currentJoint.Angle(1), Vector3d::UnitY());
        AngleAxisd rcalf(currentJoint.Angle(2), Vector3d::UnitY());
        Vector3d jointsP[4];
        jointsP[0].setZero();
        jointsP[1](0) = 0;
        jointsP[1](1) = L1 * cos(currentJoint.Angle(0));
        jointsP[1](2) = L1 * sin(currentJoint.Angle(0));
        jointsP[2](0) = jointsP[1](0) - L2 * sin(currentJoint.Angle(1));
        jointsP[2](1) = jointsP[1](1) + L2 * sin(currentJoint.Angle(0)) * cos(currentJoint.Angle(1)) + L2b * cos(currentJoint.Angle(0));
        jointsP[2](2) = jointsP[1](2) - L2 * cos(currentJoint.Angle(0)) * cos(currentJoint.Angle(1)) + L2b * sin(currentJoint.Angle(0));
        jointsP[3](0) = jointsP[2](0) - L3 * sin(currentJoint.Angle(1) + currentJoint.Angle(2));
        jointsP[3](1) = jointsP[2](1) + L3 * sin(currentJoint.Angle(0)) * cos(currentJoint.Angle(1) + currentJoint.Angle(2)) + L3b * cos(currentJoint.Angle(0));
        jointsP[3](2) = jointsP[2](2) - L3 * cos(currentJoint.Angle(0)) * cos(currentJoint.Angle(1) + currentJoint.Angle(2)) + L3b * sin(currentJoint.Angle(0));
        // 计算每个连杆的质心到单腿坐标系的位置
        this->Pcl[0] = rhip.toRotationMatrix() * Pcj[0] + jointsP[0];
        this->Pcl[1] = (rhip * rthigh).toRotationMatrix() * Pcj[1] + jointsP[1];
        this->Pcl[2] = (rhip * rthigh * rcalf).toRotationMatrix() * Pcj[2] + jointsP[2];
        this->Pcl[3] = (rhip * rthigh * rcalf).toRotationMatrix() * Pcj[3] + jointsP[3];
        // 计算每个连杆的惯性张量在单腿坐标系下的表达
        this->Icleg[0] = mI2aI(rhip.toRotationMatrix() * aI2mI(Ic[0]) * rhip.toRotationMatrix().transpose());
        this->Icleg[1] = mI2aI((rhip * rthigh).toRotationMatrix() * aI2mI(Ic[1]) * (rhip * rthigh).toRotationMatrix().transpose());
        this->Icleg[2] = mI2aI((rhip * rthigh * rcalf).toRotationMatrix() * aI2mI(Ic[2]) * (rhip * rthigh * rcalf).toRotationMatrix().transpose());
        this->Icleg[3] = mI2aI((rhip * rthigh * rcalf).toRotationMatrix() * aI2mI(Ic[3]) * (rhip * rthigh * rcalf).toRotationMatrix().transpose());
    }

    const JointS& Leg::getJointCurrent()
    {
        return currentJoint;
    }

    const LegS& Leg::getLegCurrent()
    {
        return currentLeg;
    }
}
