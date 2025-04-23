#pragma once

#include <iostream>
#include "Leg.h"
#include "kelmanFilter.h"
#include "RobotEst.h"

#define CUTOFFHZ 50

namespace Quadruped
{
    enum leg_id
    {
        LF = 0,
        RF = 1,
        LB = 2,
        RB = 3,
    };

    // 世界坐标系下整机状态，坐标系表示为{s}
    typedef struct _worldFrame
    {
        LegS leg_s[4];
        Vector4i contactEst = Vector4i::Ones(); // 接触预测
        Vector3d dist = Vector3d::Zero();       // 机身在世界坐标系下的位移
        Vector3d angVel_xyz = Vector3d::Zero(); // 角速度
        Vector3d linVel_xyz = Vector3d::Zero(); // 线速度
         Vector3d angAcc_xyz = Vector3d::Zero(); // 角加速度
        Vector3d linAcc_xyz = Vector3d::Zero(); // 线性加速度
    } worldFrame;
    // 机器人坐标系下的状态，坐标系表示为{b}
    typedef struct _bodyFrame
    {
        LegS leg_b[4];
        Vector3d Ang_xyz = Vector3d::Zero();    // 角度
        Quaterniond Quat;       //四元数
        Vector3d angVel_xyz = Vector3d::Zero(); // 角速度
        Vector3d linVel_xyz = Vector3d::Zero(); // 线速度
        Vector3d angAcc_xyz = Vector3d::Zero(); // 角加速度
        Vector3d linAcc_xyz = Vector3d::Zero(); // 线性加速度
    } bodyFrame;

    class Body
    {
    public:
        worldFrame targetWorldState;
        worldFrame currentWorldState;
        bodyFrame targetBodyState;
        bodyFrame currentBodyState;

        EstBase* est;
        Leg* legs[4];
        Vector3d leg2body[4];                  // 第一个值代表左前腿向量
        Matrix3d Rsb_c = Matrix3d::Identity(); // (当前)世界坐标系到机身坐标系的旋转矩阵映射
        Matrix3d Rsbh_c = Matrix3d::Identity();// (当前)世界坐标系到机身水平坐标系的旋转矩阵映射
        Matrix4d Tsb_c = Matrix4d::Identity(); // (当前)世界坐标系到机身坐标系的齐次变换映射
        Matrix4d Tsbh_c = Matrix4d::Identity(); // (当前)世界坐标系到机身水平坐标系的齐次变换映射
        Matrix3d Rsb_t = Matrix3d::Identity(); // (目标)世界坐标系到机身坐标系的旋转矩阵映射
        Matrix4d Tsb_t = Matrix4d::Identity(); // (目标)世界坐标系到机身坐标系的齐次变换映射
        Matrix3d dEuler2W = Matrix3d::Identity(); // 轴角速度转到欧拉角角速度
        Matrix3d Rsbf_c[4] = { Matrix3d::Identity() };// (当前)世界坐标系到足端坐标系的旋转矩阵映射
        Matrix4d Tsbf_c[4] = { Matrix4d::Identity() };// (当前)世界坐标系到足端坐标系的旋转矩阵映射
        Matrix3d Rsbw_c[4] = { Matrix3d::Identity() };// (当前)世界坐标系到轮坐标系的旋转矩阵映射
        Matrix4d Tsbw_c[4] = { Matrix4d::Identity() };// (当前)世界坐标系到轮坐标系的旋转矩阵映射

        double Mb[3] = {};      // 机器人机体质量(包括头部，中段和尾部)
        Vector<double, 6> Ib[3] = {};    // 机器人单刚体动力学机身转动惯量
        Matrix3d Mbm = Matrix3d::Zero();
        Matrix3d Ibm = Matrix3d::Zero();
        Vector3d Pbm = Vector3d::Zero();
        Vector3d Pb[3] = {};    // 机器人重心在机身坐标系下的位置
        Matrix3d M = Matrix3d::Zero(); // 机器人总质量
        Matrix3d I = Matrix3d::Zero(); // 机器人总惯量
        Vector3d P = Vector3d::Zero(); // 机器人总质心位置
        Vector3d g = Vector3d(0, 0, -9.81); // 直接初始化重力加速度向量
        double dt;                          // 控制周期

        Eigen::Matrix<double,3,4> initLegsXYPosition; // 指定右后腿初始（平面）位置

        Vector4i mixContact = Vector4i::Ones();

        // 将向量转换成角对称矩阵
        Matrix3d v3_to_m3(Vector3d _v);


    public:
        /* 构造函数，绑定四个腿部状态描述类型 */
        Body(EstBase* _est, Leg* _legObj[4], double _dt);
        /* 初始化函数，用于初始化物理参数 */
        void initParams(Vector3d _leg2body, Eigen::Matrix<double, 3, 4> _initLegsXYPosition, double _Mb[3], Vector<double, 6> _Ib[3], Vector3d _Pb[3]);
        /* 更新中性立足点 */
        void updateLegsXYPosition(Eigen::Matrix<double, 3, 4> _initLegsXYPosition);
        /* 计算齐次变换矩阵(direction为1，则是当前；为 - 1，则是目标) */
        void calTbs(int8_t direction, const Eigen::Vector3d& _slopeN);
        /* 四足运动学：改变四条腿足端的位置从而改变机器人机身的位置和姿态
           计算单腿基坐标系与机身坐标系下的足端位置转换(direction为1，则是(当前)腿->身；为-1，则是(目标)身->腿) */
        void legAndBodyPosition(int8_t direction);
        /* 计算机身坐标系与世界坐标系下的足端位置转换(世界坐标系定义为初始状态下右后腿足端位置)（direction为1，则是(当前)身->世；为 - 1，则是(目标)世->身）*/
        void bodyAndWorldFramePosition(int8_t direction);
        /* 计算足端在世界坐标系中的足端相对于机身的速度 */
        void legVelocityInWorldFrame();
        /* 计算足端再世界坐标系中的加速度 */
        void legAccInWorldFrame();
        /* 更新整机等效重心和惯量参数 */
        void updateEqBody();
        /* 接触预测 */
        void estimateContact(const Vector4i& _contact_t, const double& _t);
        void estimateContact(const Vector4i& _contact_t, const double& _t, const Matrix3d& _slope);

        /* 更新目标位姿 */
        void updateBodyTargetPos(Vector3d _angle, Vector3d _position);
        /* 更新惯导姿态 */
        void updateBodyImu(Vector3d _imuRPY);
        void updateBodyImu(Quaterniond _imyQ);
        /* 更新陀螺仪角速度 */
        void updateBodyGyro(Vector3d _gyro);
        void updateBodyGyroAcc(Vector3d _gyroAcc);
        /* 更新加速度计加速度 */
        void updateBodyAcc(Vector3d _acc);
        /* 更新目标四足点(地面) */
        void updateTargetFootPoint(Eigen::Matrix<double,3,4> _point);
        /* 更新目标四足速度 */
        void updateTargetFootVel(Eigen::Matrix<double, 3, 4> _vel);

        /* 获取运动学四足点位 */
        Eigen::Matrix<double, 3, 4> getFKFeetPos();
        Eigen::Vector3d getFKFeetPos(int id);
        Eigen::Matrix<double, 3, 4> getFKFeetPosB();
        Eigen::Vector3d getFKFeetPosB(int id);
        Eigen::Matrix<double, 3, 4> getFKFeetVel();
        Eigen::Vector3d getFKFeetVel(int i);

        /* 一些数学计算函数 */
        Eigen::Matrix3d quatToRot(const Vector4d& _quat);
        Eigen::Vector3d rotMatToEulerZYX(const Matrix3d& R);
        Eigen::Vector<double, 6> parallelAxis(const Vector3d& dstAxis, const Vector3d& oriAxis, const Vector<double, 6>& oriIm, const double oriM);// 平行轴定理
        Eigen::Matrix3d aI2mI(const Vector<double, 6>& oriIm);
        Eigen::Matrix3d eulerVelToRotVel(const Vector3d& _euler);
        Eigen::Matrix3d axisToRotationMatrixQ(const Eigen::Vector3d& _u, const Eigen::Vector3d& _v);
        Eigen::Matrix3d axisToRotationMatrixR(const Eigen::Vector3d& _u, const Eigen::Vector3d& _v);
    };

    Body::Body(EstBase* _est,Leg* _legObj[4], double _dt)
    {
        est = _est;
        for (int i = 0; i < 4; i++)
        {
            this->legs[i] = _legObj[i];
        }
        this->dt = _dt;
    }

    void Body::initParams(Vector3d _leg2body, Eigen::Matrix<double,3,4> _initLegsXYPosition, double _Mb[3], Vector<double, 6> _Ib[3], Vector3d _Pb[3])
    {
        // 腿部坐标系到机身坐标系的位置偏移
        this->leg2body[LF] = _leg2body;
        Vector3d tmp = _leg2body;
        tmp(1) = tmp(1) * -1;
        this->leg2body[RF] = tmp;
        tmp(0) = tmp(0) * -1;
        this->leg2body[RB] = tmp;
        tmp(1) = tmp(1) * -1;
        this->leg2body[LB] = tmp;
        // 对于对称的步态定义一个初始参数，描述的是中性立足点
        this->initLegsXYPosition = _initLegsXYPosition;
      /*  this->M = Eigen::Map<Matrix3d>(_M.data(),3,3);
        this->Ib = Eigen::Map<Matrix3d>(_Ib.data(),3,3);*/
        for (int i = 0; i < 3; i++)
        {
            this->Mb[i] = _Mb[i];
            this->Ib[i] = _Ib[i];
            this->Pb[i] = _Pb[i];
        }
    }

    void Body::updateLegsXYPosition(Eigen::Matrix<double, 3, 4> _initLegsXYPosition)
    {
        this->initLegsXYPosition = _initLegsXYPosition;
    }

    Matrix3d Body::v3_to_m3(Vector3d _v)
    {
        Matrix3d m = Matrix3d::Zero();
        m(0, 1) = -_v(2);
        m(0, 2) = _v(1);
        m(1, 0) = _v(2);
        m(1, 2) = -_v(0);
        m(2, 0) = -_v(1);
        m(2, 1) = _v(0);
        return m;
    }

    void Body::calTbs(int8_t direction, const Eigen::Vector3d& _slopeN)
    {
        if (direction == 1)
        {
            // 使用四元数得到变换矩阵
            /*Eigen::Matrix3d rotation_c = quatToRot(currentBodyState.Quat);*/
            Eigen::Matrix3d rotation_c = currentBodyState.Quat.toRotationMatrix();
            Rsb_c = rotation_c;
            currentBodyState.Ang_xyz = rotMatToEulerZYX(Rsb_c);
            Eigen::AngleAxisd rotationz_c(currentBodyState.Ang_xyz(2), Eigen::Vector3d::UnitZ());
            Eigen::AngleAxisd rotationy_c(currentBodyState.Ang_xyz(1), Eigen::Vector3d::UnitY());
            Eigen::AngleAxisd rotationx_c(currentBodyState.Ang_xyz(0), Eigen::Vector3d::UnitX());
            Rsbh_c = rotationz_c.toRotationMatrix();
            Tsb_c.block<3, 3>(0, 0) = Rsb_c;
            Tsb_c.block<3, 1>(0, 3) = currentWorldState.dist;
            Tsbh_c.block<3, 3>(0, 0) = Rsbh_c;
            Tsbh_c.block<3, 1>(0, 3) = currentWorldState.dist;
            this->dEuler2W = this->eulerVelToRotVel(currentBodyState.Ang_xyz);
            // 更新其他世界坐标系下的变量
            currentWorldState.angVel_xyz = Rsb_c * currentBodyState.angVel_xyz;
            currentWorldState.linAcc_xyz = Rsb_c * currentBodyState.linAcc_xyz;
            est->updateTsb(Tsb_c);
            // 足端姿态正运动学
            for (int i = 0; i < 4; i++)
            {
                Eigen::AngleAxisd rotation_hip(this->legs[i]->currentJoint.Angle(0), Eigen::Vector3d::UnitX());
                Eigen::AngleAxisd rotation_thigh(this->legs[i]->currentJoint.Angle(1), Eigen::Vector3d::UnitY());
                Eigen::AngleAxisd rotation_calf(this->legs[i]->currentJoint.Angle(2), Eigen::Vector3d::UnitY());
                this->Rsbf_c[i] = this->Rsb_c * (rotation_hip * rotation_thigh * rotation_calf).toRotationMatrix();
                this->Tsbf_c[i].block(0, 0, 3, 3) = this->Rsbf_c[i];
                this->Tsbf_c[i].block(0, 3, 3, 1) = currentWorldState.dist + this->currentWorldState.leg_s[i].Position;
            }
            // 轮坐标系简化计算，直接通过地形得到z轴，足端姿态正运动学得到y轴，叉乘得到x轴，然后用等效轴角计算旋转矩阵
            for (int i = 0; i < 4; i++)
            {
                Eigen::Vector3d zAxis = _slopeN;
                Eigen::Vector3d yAxis = this->Rsbf_c[i] * Eigen::Vector3d(0, 1., 0);
                Eigen::Vector3d xAxis = yAxis.cross(zAxis);
                xAxis.normalize();
                Eigen::Vector3d refNormal(1., 0, 0);
                this->Rsbw_c[i] = this->axisToRotationMatrixQ(refNormal, xAxis);
                this->Tsbw_c[i].block(0, 0, 3, 3) = this->Rsbw_c[i];
                this->Tsbw_c[i].block(0, 3, 3, 1) = currentWorldState.dist + this->currentWorldState.leg_s[i].Position;
                //std::cout << i << ": " << this->Rsbw_c[i].transpose() << std::endl;
                /*std::cout << xAxis << std::endl;
                std::cout << _slopeN << std::endl;*/
            }
        }
        else if (direction == -1)
        {
            // 目标
            // 三轴欧拉角旋转矩阵
            Eigen::AngleAxisd rotationx_t(targetBodyState.Ang_xyz(0), Eigen::Vector3d::UnitX());
            Eigen::AngleAxisd rotationy_t(targetBodyState.Ang_xyz(1), Eigen::Vector3d::UnitY());
            Eigen::AngleAxisd rotationz_t(targetBodyState.Ang_xyz(2), Eigen::Vector3d::UnitZ());
            Eigen::Matrix3d rotation_t = (rotationz_t * rotationy_t * rotationx_t).toRotationMatrix();
            // 构建齐次变换矩阵
            Rsb_t = rotation_t;
            Tsb_t.block<3, 3>(0, 0) = rotation_t;
            Tsb_t.block<3, 1>(0, 3) = targetWorldState.dist;
            // 更新其他世界坐标系下的变量
            // targetBodyState.angVel_xyz = rotation_t.inverse() * targetWorldState.angVel_xyz;
            /*targetBodyState.linVel_xyz = rotation_t.inverse() * targetWorldState.linVel_xyz;*/
            // targetBodyState.angAcc_xyz = rotation_t.inverse() * targetWorldState.angAcc_xyz;
            /*targetBodyState.linAcc_xyz = rotation_t.inverse() * targetWorldState.linAcc_xyz;*/
        }
    }

    void Body::legAndBodyPosition(int8_t direction)
    {
        if (direction == 1)
        {
            currentBodyState.leg_b[LF].Position = legs[LF]->currentLeg.Position + this->leg2body[LF];
            currentBodyState.leg_b[RF].Position = legs[RF]->currentLeg.Position + this->leg2body[RF];
            currentBodyState.leg_b[RB].Position = legs[RB]->currentLeg.Position + this->leg2body[RB];
            currentBodyState.leg_b[LB].Position = legs[LB]->currentLeg.Position + this->leg2body[LB];
            currentBodyState.leg_b[LF].Force = legs[LF]->currentLeg.Force;
            currentBodyState.leg_b[RF].Force = legs[RF]->currentLeg.Force;
            currentBodyState.leg_b[RB].Force = legs[RB]->currentLeg.Force;
            currentBodyState.leg_b[LB].Force = legs[LB]->currentLeg.Force;
        }
        else if (direction == -1)
        {
            legs[LF]->targetLeg.Position = targetBodyState.leg_b[LF].Position - this->leg2body[LF];
            legs[RF]->targetLeg.Position = targetBodyState.leg_b[RF].Position - this->leg2body[RF];
            legs[RB]->targetLeg.Position = targetBodyState.leg_b[RB].Position - this->leg2body[RB];
            legs[LB]->targetLeg.Position = targetBodyState.leg_b[LB].Position - this->leg2body[LB];
        }
    }

    void Body::bodyAndWorldFramePosition(int8_t direction)
    {
        for (int i = 0; i < 4; i++)
        {
            if (direction == 1)
            {
                Vector4d Pbi(0, 0, 0, 1);
                Pbi.block<3, 1>(0, 0) = currentBodyState.leg_b[i].Position;
                Vector4d Psi = this->Tsb_c * Pbi;
                currentWorldState.leg_s[i].Position = Psi.block<3, 1>(0, 0);
                currentWorldState.leg_s[i].Force = this->Rsb_c * currentBodyState.leg_b[i].Force;
            }
            else if (direction == -1)
            {
                Vector4d Psi(0, 0, 0, 1);
                Psi.block<3, 1>(0, 0) = targetWorldState.leg_s[i].Position;
                Vector4d Pbi = this->Tsb_t.inverse() * Psi;
                targetBodyState.leg_b[i].Position = Pbi.block<3, 1>(0, 0);
            }
        }
    }

    void Body::legVelocityInWorldFrame()
    {
        // 将机身旋转角速度向量转换成反对称矩阵
        Matrix3d w = v3_to_m3(currentBodyState.angVel_xyz);
        for (int i = 0; i < 4; i++)
        {
            // 更新足端相对于机身坐标系的速度
            currentBodyState.leg_b[i].Velocity = legs[i]->currentLeg.Velocity;
            currentBodyState.leg_b[i].VelocityW = legs[i]->currentLeg.VelocityW;
            currentBodyState.leg_b[i].VelocityG = legs[i]->currentLeg.VelocityG;
            // 计算足端相对于世界坐标系的速度
            currentWorldState.leg_s[i].Velocity = Rsb_c * (w * currentBodyState.leg_b[i].Position + legs[i]->currentLeg.Velocity);
            currentWorldState.leg_s[i].VelocityW = Rsbh_c * legs[i]->currentLeg.VelocityW;
            currentWorldState.leg_s[i].VelocityG = currentWorldState.leg_s[i].Velocity + currentWorldState.leg_s[i].VelocityW;
        }
    }

    void Body::legAccInWorldFrame()
    {
        /*static LPF_SecondOrder_Classdef legAccF[4][3] = { {LPF_SecondOrder_Classdef(CUTOFFHZ,500),LPF_SecondOrder_Classdef(CUTOFFHZ,500),LPF_SecondOrder_Classdef(CUTOFFHZ,500)},\
                                                            {LPF_SecondOrder_Classdef(CUTOFFHZ, 500), LPF_SecondOrder_Classdef(CUTOFFHZ, 500), LPF_SecondOrder_Classdef(CUTOFFHZ, 500)},\
                                                            {LPF_SecondOrder_Classdef(CUTOFFHZ, 500), LPF_SecondOrder_Classdef(CUTOFFHZ, 500), LPF_SecondOrder_Classdef(CUTOFFHZ, 500)},\
                                                            {LPF_SecondOrder_Classdef(CUTOFFHZ, 500), LPF_SecondOrder_Classdef(CUTOFFHZ, 500), LPF_SecondOrder_Classdef(CUTOFFHZ, 500)} };*/
        // 将机身旋转角速度向量转换成反对称矩阵
        Matrix3d w = v3_to_m3(currentBodyState.angVel_xyz);
        Matrix3d w_dot = v3_to_m3(currentBodyState.angAcc_xyz);
        for (int i = 0; i < 4; i++)
        {
            // 更新足端相对于机身坐标系的加速度
            currentBodyState.leg_b[i].Acc = legs[i]->currentLeg.Acc;
            // 更新足端相对于世界坐标系的加速度
            Vector3d Acc = currentWorldState.linAcc_xyz + Rsb_c * (w_dot * currentBodyState.leg_b[i].Position + w * w * currentBodyState.leg_b[i].Position + legs[i]->currentLeg.Acc);
            currentWorldState.leg_s[i].Acc = Acc;
            /*currentWorldState.leg_s[i].Acc(0) = legAccF[i][0].f(Acc(0));
            currentWorldState.leg_s[i].Acc(1) = legAccF[i][1].f(Acc(1));
            currentWorldState.leg_s[i].Acc(2) = legAccF[i][2].f(Acc(2));*/
        }
    }

    void Body::updateEqBody()
    {
        // 首先计算整机重量惯量参数
        double _M = 0;
        Vector3d _P = Vector3d::Zero();
        Vector<double, 6> _I = Vector<double, 6>::Zero();
        for (int i = 0; i < 3; i++)
        {
            _M += this->Mb[i];
            _P += this->Mb[i] * this->Pb[i];
        }

        /* 只考虑机身单刚体 */
        this->Pbm = _P / _M;// 只包含机身的的重心位置
        this->Mbm = _M * Matrix3d::Identity();// 只包含机身的重量
        // 只包含机身的惯量
        Vector<double, 6> _Ibm = Vector<double, 6>::Zero();
        for (int i = 0; i < 3; i++)
        {
            _Ibm += parallelAxis(this->Pbm, Pb[i], Ib[i], Mb[i]);
        }
        this->Ibm = aI2mI(_Ibm);

        /* 考虑全身 */
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                _M += this->legs[i]->Mc[j];
                _P += this->legs[i]->Mc[j] * (this->legs[i]->Pcl[j] + this->leg2body[i]);
            }
        }
        _P = _P / _M;// 计算整机重心位置
        for (int i = 0; i < 3; i++)
        {
            _I += parallelAxis(_P, Pb[i], Ib[i], Mb[i]);
        }
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                _I += parallelAxis(_P, this->legs[i]->Pcl[j] + this->leg2body[i], this->legs[i]->Icleg[j], this->legs[i]->Mc[j]);
            }
        }
        // 得到整机刚体参数
        this->P = _P;
        this->M = _M * Matrix3d::Identity();
        this->I = aI2mI(_I);

    }

    void Body::estimateContact(const Vector4i& _contact_t, const double& _t)
    {
        static LPF_SecondOrder_Classdef extForceF[4][3] = { {LPF_SecondOrder_Classdef(CUTOFFHZ,500),LPF_SecondOrder_Classdef(CUTOFFHZ,500),LPF_SecondOrder_Classdef(CUTOFFHZ,500)},\
                                                            {LPF_SecondOrder_Classdef(CUTOFFHZ, 500), LPF_SecondOrder_Classdef(CUTOFFHZ, 500), LPF_SecondOrder_Classdef(CUTOFFHZ, 500)},\
                                                            {LPF_SecondOrder_Classdef(CUTOFFHZ, 500), LPF_SecondOrder_Classdef(CUTOFFHZ, 500), LPF_SecondOrder_Classdef(CUTOFFHZ, 500)},\
                                                            {LPF_SecondOrder_Classdef(CUTOFFHZ, 500), LPF_SecondOrder_Classdef(CUTOFFHZ, 500), LPF_SecondOrder_Classdef(CUTOFFHZ, 500)} };
        static int vacantCheck[4] = { 1,1,1,1 };
        static double last_t[4] = { 0 };
        static Vector4i lastContact_c = Vector4i::Ones();
        static Vector4i lastContact_t = Vector4i::Ones();

        this->targetWorldState.contactEst = _contact_t;
        for (int i = 0; i < 4; i++)
        {
            Vector3d extForce = this->legs[i]->Mc[3] * (this->currentWorldState.leg_s[i].Acc/* - Vector3d(0, 0, -9.81)*/) - this->currentWorldState.leg_s[i].Force;
            //this->currentWorldState.leg_s[i].extForce = extForce;
            this->currentWorldState.leg_s[i].extForce(0) = extForceF[i][0].f(extForce(0));
            this->currentWorldState.leg_s[i].extForce(1) = extForceF[i][1].f(extForce(1));
            this->currentWorldState.leg_s[i].extForce(2) = extForceF[i][2].f(extForce(2));
            if (vacantCheck[i] == 1 && this->currentWorldState.leg_s[i].extForce(2) < 100.)
            {
                this->currentWorldState.contactEst(i) = 0;
                vacantCheck[i] = 0;
            }
            if (this->currentWorldState.leg_s[i].Acc(2) > 225.)
            {
                this->currentWorldState.contactEst(i) = 1;
                last_t[i] = _t;// 记录接触的时刻
            }
            if (this->currentWorldState.contactEst(i) == 0 && this->currentWorldState.leg_s[i].extForce(2) > 125.)
            {
                this->currentWorldState.contactEst(i) = 1;
            }
            if (vacantCheck[i] == 0 && this->currentWorldState.contactEst(i) == 1 && (_t - last_t[i]) > 0.1)
            {
                vacantCheck[i] = 1;
            }

            // 混合接触相位，上升沿参考被动检测，下降沿参考步态周期规划
            if (mixContact(i) == 0)
            {
                if (this->currentWorldState.contactEst(i) == 1 && lastContact_c(i) == 0)
                {
                    mixContact(i) = 1;
                }
                if (this->targetWorldState.contactEst(i) == 1)
                {
                    mixContact(i) = 1;
                }
            }
            else
            {
                if (this->targetWorldState.contactEst(i) == 0 && lastContact_t(i) == 1)
                {
                    mixContact(i) = 0;
                }
            }
        }
        lastContact_c = this->currentWorldState.contactEst;
        lastContact_t = this->targetWorldState.contactEst;
    }

    void Body::estimateContact(const Vector4i& _contact_t, const double& _t, const Matrix3d& _slope)
    {
        static LPF_SecondOrder_Classdef extForceF[4][3] = { {LPF_SecondOrder_Classdef(CUTOFFHZ,500),LPF_SecondOrder_Classdef(CUTOFFHZ,500),LPF_SecondOrder_Classdef(CUTOFFHZ,500)},\
                                                            {LPF_SecondOrder_Classdef(CUTOFFHZ, 500), LPF_SecondOrder_Classdef(CUTOFFHZ, 500), LPF_SecondOrder_Classdef(CUTOFFHZ, 500)},\
                                                            {LPF_SecondOrder_Classdef(CUTOFFHZ, 500), LPF_SecondOrder_Classdef(CUTOFFHZ, 500), LPF_SecondOrder_Classdef(CUTOFFHZ, 500)},\
                                                            {LPF_SecondOrder_Classdef(CUTOFFHZ, 500), LPF_SecondOrder_Classdef(CUTOFFHZ, 500), LPF_SecondOrder_Classdef(CUTOFFHZ, 500)} };
        static int vacantCheck[4] = { 1,1,1,1 };
        static double last_t[4] = { 0 };
        static Vector4i lastContact_c = Vector4i::Ones();
        static Vector4i lastContact_t = Vector4i::Ones();

        this->targetWorldState.contactEst = _contact_t;
        for (int i = 0; i < 4; i++)
        {
            Vector3d extForce = this->legs[i]->Mc[3] * (this->currentWorldState.leg_s[i].Acc/* - Vector3d(0, 0, -9.81)*/) - this->currentWorldState.leg_s[i].Force;
            //this->currentWorldState.leg_s[i].extForce = extForce;
            this->currentWorldState.leg_s[i].extForce(0) = extForceF[i][0].f(extForce(0));
            this->currentWorldState.leg_s[i].extForce(1) = extForceF[i][1].f(extForce(1));
            this->currentWorldState.leg_s[i].extForce(2) = extForceF[i][2].f(extForce(2));
            if (vacantCheck[i] == 1 && (_slope.transpose() * this->currentWorldState.leg_s[i].extForce)(2) < 100.)
            {
                this->currentWorldState.contactEst(i) = 0;
                vacantCheck[i] = 0;
            }
            if ((_slope.transpose() * this->currentWorldState.leg_s[i].Acc)(2) > 225.)
            {
                this->currentWorldState.contactEst(i) = 1;
                last_t[i] = _t;// 记录接触的时刻
            }
            if (this->currentWorldState.contactEst(i) == 0 && (_slope.transpose() * this->currentWorldState.leg_s[i].extForce)(2) > 125.)
            {
                this->currentWorldState.contactEst(i) = 1;
            }
            if (vacantCheck[i] == 0 && this->currentWorldState.contactEst(i) == 1 && (_t - last_t[i]) > 0.1)
            {
                vacantCheck[i] = 1;
            }

            // 混合接触相位，上升沿参考被动检测，下降沿参考步态周期规划
            if (mixContact(i) == 0)
            {
                if (this->currentWorldState.contactEst(i) == 1 && lastContact_c(i) == 0)
                {
                    mixContact(i) = 1;
                }
                if (this->targetWorldState.contactEst(i) == 1)
                {
                    mixContact(i) = 1;
                }
            }
            else
            {
                if (this->targetWorldState.contactEst(i) == 0 && lastContact_t(i) == 1)
                {
                    mixContact(i) = 0;
                }
            }
        }
        lastContact_c = this->currentWorldState.contactEst;
        lastContact_t = this->targetWorldState.contactEst;
    }

    void Body::updateBodyTargetPos(Vector3d _angle, Vector3d _position)
    {
        targetBodyState.Ang_xyz = _angle;
        targetWorldState.dist = _position;
    }

    void Body::updateBodyImu(Vector3d _imuRPY)
    {
        currentBodyState.Ang_xyz = _imuRPY;
    }

    void Body::updateBodyImu(Quaterniond _imuQ)
    {
        currentBodyState.Quat = _imuQ;
    }

    void Body::updateBodyGyro(Vector3d _gyro)
    {
        currentBodyState.angVel_xyz = _gyro;
    }

    void Body::updateBodyGyroAcc(Vector3d _gyroAcc)
    {
        currentBodyState.angAcc_xyz = _gyroAcc;
    }

    void Body::updateBodyAcc(Vector3d _acc)
    {
        currentBodyState.linAcc_xyz = _acc;
    }

    void Body::updateTargetFootPoint(Eigen::Matrix<double, 3, 4> _point)
    {
        for (int i = 0; i < 4; i++)
        {
            targetWorldState.leg_s[i].Position = _point.col(i);
            Vector4d Psi(0, 0, 0, 1);
            Psi.block<3, 1>(0, 0) = targetWorldState.leg_s[i].Position;
            Vector4d Pbi = this->Tsb_c.inverse() * Psi;
            targetBodyState.leg_b[i].Position = Pbi.block<3, 1>(0, 0);
        }
    }

    void Body::updateTargetFootVel(Eigen::Matrix<double, 3, 4> _vel)
    {
        for (int i = 0; i < 4; i++)
        {
            //此处是指足端相对于机身的速度在世界坐标系中的表达，因此需要减去整机速度
            targetWorldState.leg_s[i].Velocity = _vel.col(i) - est->getEstBodyVelS();
            targetBodyState.leg_b[i].Velocity = Rsb_c.transpose() * _vel.col(i);
            legs[i]->setTargetLegVelocity(targetBodyState.leg_b[i].Velocity);
        }
    }

    Eigen::Matrix<double, 3, 4> Body::getFKFeetPos()
    {
        Eigen::Matrix<double, 3, 4> out;
        for (int i = 0; i < 4; i++)
        {
            out.block<3, 1>(0, i) = currentWorldState.leg_s[i].Position;
        }
        return out;
    }

    Eigen::Vector3d Body::getFKFeetPos(int id)
    {
        Eigen::Vector3d out;
        out = currentWorldState.leg_s[id].Position;
        return out;
    }

    Eigen::Matrix<double, 3, 4> Body::getFKFeetPosB()
    {
        Eigen::Matrix<double, 3, 4> out;
        for (int i = 0; i < 4; i++)
        {
            out.block<3, 1>(0, i) = currentBodyState.leg_b[i].Position;
        }
        return out;
    }

    Eigen::Vector3d Body::getFKFeetPosB(int id)
    {
        Eigen::Vector3d out;
        out = currentBodyState.leg_b[id].Position;
        return out;
    }

    Eigen::Matrix<double, 3, 4> Body::getFKFeetVel()
    {
        Eigen::Matrix<double, 3, 4> out;
        for (int i = 0; i < 4; i++)
        {
            out.block<3, 1>(0, i) = currentWorldState.leg_s[i].VelocityG;
        }
        return out;
    }

    Eigen::Vector3d Body::getFKFeetVel(int id)
    {
        Eigen::Vector3d out;
        out = currentWorldState.leg_s[id].VelocityG;
        return out;
    }

    Eigen::Matrix3d Body::quatToRot(const Vector4d& _quat)
    {
        double e0 = _quat(0);
        double e1 = _quat(1);
        double e2 = _quat(2);
        double e3 = _quat(3);

        Matrix3d R;
        R << 1 - 2 * (e2 * e2 + e3 * e3), 2 * (e1 * e2 - e0 * e3),
            2 * (e1 * e3 + e0 * e2), 2 * (e1 * e2 + e0 * e3),
            1 - 2 * (e1 * e1 + e3 * e3), 2 * (e2 * e3 - e0 * e1),
            2 * (e1 * e3 - e0 * e2), 2 * (e2 * e3 + e0 * e1),
            1 - 2 * (e1 * e1 + e2 * e2);
        return R;
    }

    Eigen::Vector3d Body::rotMatToEulerZYX(const Matrix3d& R) {
        Vector3d rpy;// 顺序为ZYX
        constexpr double PI = 3.14159265358979323846;
        constexpr double EPSILON = 1e-6;
        // 检查万向节锁（俯仰角为±90度）
        if (std::abs(R(2,0)) > 1.0 - EPSILON) {
            rpy(1) = (R(2,0) > 0) ? PI / 2 : -PI / 2;
            rpy(0) = 0.0; // 固定横滚角为0

            // 计算偏航角（通过矩阵不同元素组合）
            rpy(2) = std::atan2((R(2,0) > 0) ? R(1,2) : -R(1,2),R(1,1));
        }
        else {
            rpy(2) = std::atan2(R(1, 0), R(0, 0));
            rpy(1) = std::asin(-R(2, 0));
            rpy(0) = std::atan2(R(2, 1), R(2, 2));
        }
        return rpy;
    }

    Eigen::Vector<double, 6> Body::parallelAxis(const Vector3d& dstAxis, const Vector3d& oriAxis, const Vector<double, 6>& oriIm, const double oriM)
    {
        Vector<double, 6> dstIm;
        Vector3d eAxis = oriAxis - dstAxis;
        dstIm(0) = oriIm(0) + oriM * (pow(eAxis(1), 2) + pow(eAxis(2), 2));
        dstIm(1) = oriIm(1) + oriM * (pow(eAxis(0), 2) + pow(eAxis(2), 2));
        dstIm(2) = oriIm(2) + oriM * (pow(eAxis(0), 2) + pow(eAxis(1), 2));
        dstIm(3) = oriIm(3) - oriM * eAxis(0) * eAxis(1);
        dstIm(4) = oriIm(4) - oriM * eAxis(0) * eAxis(2);
        dstIm(5) = oriIm(5) - oriM * eAxis(1) * eAxis(2);
        return dstIm;
    }

    Eigen::Matrix3d Body::aI2mI(const Vector<double, 6>& oriIa)
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

    Eigen::Matrix3d Body::eulerVelToRotVel(const Vector3d& _euler)
    {
        Matrix3d t = Matrix3d::Zero();
        t(0, 0) = 1;
        t(0, 2) = -std::sin(_euler(1));
        t(1, 1) = std::cos(_euler(0));
        t(1, 2) = std::cos(_euler(1)) * std::sin(_euler(0));
        t(2, 1) = -std::sin(_euler(0));
        t(2, 2) = std::cos(_euler(1)) * std::cos(_euler(0));
        return t;
    }

    Eigen::Matrix3d Body::axisToRotationMatrixQ(const Eigen::Vector3d& _u, const Eigen::Vector3d& _v) {
        // 归一化输入向量
        Eigen::Vector3d u_norm = _u.normalized();
        Eigen::Vector3d v_norm = _v.normalized();

        // 计算旋转轴
        Eigen::Vector3d k = u_norm.cross(v_norm);
        double k_norm = k.norm();
        double dot = u_norm.dot(v_norm);

        // 处理共线情况（旋转角为 0 或 180°）
        if (k_norm < 1e-3) {
            // 如果 u 和 v 方向相同，返回单位矩阵
            if (dot > 0) {
                return Eigen::Matrix3d::Identity();
            }
            // 如果 u 和 v 方向相反，返回绕任意垂直轴的 180° 旋转
            else {
                // 找到一个与 u 垂直的向量
                //Eigen::Vector3d perpendicular = u_norm.unitOrthogonal();
                Eigen::Vector3d perpendicular(0, 0, 1.);
                Eigen::AngleAxisd rotation(M_PI, perpendicular);
                //std::cout << "stuck: " << std::endl;
                return rotation.toRotationMatrix();
            }
        }
        // 不能接近pi，否则已经出现问题
        if (k_norm < 0.15)
        {
            if (dot < 0)
            {
                Eigen::Vector3d perpendicular(0, 0, 1.);
                Eigen::AngleAxisd rotation(M_PI, perpendicular);
                //std::cout << "stuck: " << std::endl;
                return rotation.toRotationMatrix();
            }
        }

        // 归一化旋转轴
        k.normalize();

        // 计算旋转角
        double cos_theta = u_norm.dot(v_norm);
        double theta = acos(cos_theta);

        // 构造四元数
        Eigen::Quaterniond q;
        q = Eigen::AngleAxisd(theta, k);

        // 转换为旋转矩阵
        return q.toRotationMatrix();
    }

    Eigen::Matrix3d Body::axisToRotationMatrixR(const Eigen::Vector3d& _u, const Eigen::Vector3d& _v)
    {
        // 归一化输入向量
        Eigen::Vector3d u_norm = _u.normalized();
        Eigen::Vector3d v_norm = _v.normalized();

        // 计算旋转轴
        Eigen::Vector3d k = u_norm.cross(v_norm);
        double k_norm = k.norm();
        double dot = u_norm.dot(v_norm);

        // 处理共线情况（旋转角为 0 或 180°）
        if (k_norm < 1e-3) {
            // 如果 u 和 v 方向相同，返回单位矩阵
            if (dot > 0) {
                return Eigen::Matrix3d::Identity();
            }
            // 如果 u 和 v 方向相反，返回绕任意垂直轴的 180° 旋转
            else {
                // 找到一个与 u 垂直的向量
                //Eigen::Vector3d perpendicular = u_norm.unitOrthogonal();
                Eigen::Vector3d perpendicular(0, 0, 1.);
                Eigen::AngleAxisd rotation(M_PI, perpendicular);
                //std::cout << "stuck: " << std::endl;
                return rotation.toRotationMatrix();
            }
        }

        if (k_norm < 0.15)
        {
            if (dot < 0)
            {
                Eigen::Vector3d perpendicular(0, 0, 1.);
                Eigen::AngleAxisd rotation(M_PI, perpendicular);
                //std::cout << "stuck: " << std::endl;
                return rotation.toRotationMatrix();
            }
        }

        // 归一化旋转轴
        k.normalize();

        double costheta = dot;
        double sintheta = k_norm;
        Matrix3d K;
        K << 0, -k(2), k(1), k(2), 0, -k(0), -k(1), k(0), 0;
        Matrix3d R = Matrix3d::Identity() + sintheta * K + (1 - costheta) * K * K;
        return R;
    }
}
