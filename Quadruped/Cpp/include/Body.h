#pragma once

#include <iostream>
#include "Leg.h"
#include "kelmanFilter.h"
#include "RobotEst.h"


namespace Quadruped
{
    /*template<typename T>
    inline T windowFunc(const T x, const T windowRatio, const T xRange = 1.0, const T yRange = 1.0) {
        if ((x < 0) || (x > xRange)) {
            std::cout << "[ERROR][windowFunc] The x=" << x << ", which should between [0, xRange]" << std::endl;
        }
        if ((windowRatio <= 0) || (windowRatio >= 0.5)) {
            std::cout << "[ERROR][windowFunc] The windowRatio=" << windowRatio << ", which should between [0, 0.5]" << std::endl;
        }

        if (x / xRange < windowRatio) {
            return x * yRange / (xRange * windowRatio);
        }
        else if (x / xRange > 1 - windowRatio) {
            return yRange * (xRange - x) / (xRange * windowRatio);
        }
        else {
            return yRange;
        }
    }*/

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
        Vector3d dist = Vector3d::Zero();       // 机身在世界坐标系下的位移
        Vector3d angVel_xyz = Vector3d::Zero(); // 角速度
        Vector3d linVel_xyz = Vector3d::Zero(); // 线速度
        // Vector3d angAcc_xyz = Vector3d::Zero(); // 角加速度
        Vector3d linAcc_xyz = Vector3d::Zero(); // 线性加速度
    } worldFrame;
    // 机器人坐标系下的状态，坐标系表示为{b}
    typedef struct _bodyFrame
    {
        LegS leg_b[4];
        Vector3d Ang_xyz = Vector3d::Zero();    // 角度
        Vector4d Quat = Vector4d::Zero();       //四元数
        Vector3d angVel_xyz = Vector3d::Zero(); // 角速度
        Vector3d linVel_xyz = Vector3d::Zero(); // 线速度
        //Vector3d angAcc_xyz = Vector3d::Zero(); // 角加速度
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
        Matrix4d Tsb_c = Matrix4d::Identity(); // (当前)世界坐标系到机身坐标系的齐次变换映射
        Matrix3d Rsb_t = Matrix3d::Identity(); // (目标)世界坐标系到机身坐标系的旋转矩阵映射
        Matrix4d Tsb_t = Matrix4d::Identity(); // (目标)世界坐标系到机身坐标系的齐次变换映射

        double Mb[3] = {};      // 机器人机体质量(包括头部，中段和尾部)
        Vector<double, 6> Ib[3] = {};    // 机器人单刚体动力学机身转动惯量
        Vector3d Pb[3] = {};    // 机器人重心在机身坐标系下的位置
        Matrix3d M = Matrix3d::Zero(); // 机器人总质量
        Matrix3d I = Matrix3d::Zero(); // 机器人总惯量
        Vector3d P = Vector3d::Zero(); // 机器人总质心位置
        Vector3d g = Vector3d(0, 0, -9.81); // 直接初始化重力加速度向量
        double dt;                          // 控制周期

        Eigen::Matrix<double, 6, 12> dynamicLeft = Eigen::Matrix<double, 6, 12>::Zero();
        Eigen::Matrix<double, 6, 6> dynamicRight = Eigen::Matrix<double, 6, 6>::Zero();

        Vector3d initRbLegXYPosition; // 指定右后腿初始（平面）位置

        // 将向量转换成角对称矩阵
        Matrix3d v3_to_m3(Vector3d _v);


    public:
        // 构造函数，绑定四个腿部状态描述类型
        Body(EstBase* _est, Leg* _legObj[4], double _dt);
        // 初始化函数，用于初始化物理参数
        void initParams(Vector3d _leg2body, Vector3d _initRbLegXYPosition, double _Mb[3], Vector<double, 6> _Ib[3], Vector3d _Pb[3]);
        // 计算齐次变换矩阵(direction为1，则是当前；为-1，则是目标)
        void calTbs(int8_t direction);
        // 四足运动学：改变四条腿足端的位置从而改变机器人机身的位置和姿态
        // 计算单腿基坐标系与机身坐标系下的足端位置转换(direction为1，则是(当前)腿->身；为-1，则是(目标)身->腿)
        void legAndBodyPosition(int8_t direction);
        // 计算机身坐标系与世界坐标系下的足端位置转换(世界坐标系定义为初始状态下右后腿足端位置)（direction为1，则是(当前)身->世；为-1，则是(目标)世->身）
        void bodyAndWorldFramePosition(int8_t direction);
        // 计算足端在世界坐标系中的足端速度
        void legVelocityInWorldFrame();
        // 更新机器人动力学方程（描述为 left*[f] = right）
        void updateDynamic();

        // 更新目标位姿
        void updateBodyTargetPos(Vector3d _angle, Vector3d _position);
        // 更新惯导姿态
        void updateBodyImu(Vector3d _imuRPY);
        void updateBodyImu(Vector4d _imyQ);
        // 更新陀螺仪角速度
        void updateBodyGyro(Vector3d _gyro);
        // 更新加速度计加速度
        void updateBodyAcc(Vector3d _acc);
        // 更新目标四足点(地面)
        void updateTargetFootPoint(Eigen::Matrix<double,3,4> _point);
        // 更新目标四足速度
        void updateTargetFootVel(Eigen::Matrix<double, 3, 4> _vel);

        // 获取运动学四足点位
        Eigen::Matrix<double, 3, 4> getFKFeetPos();
        Eigen::Vector3d getFKFeetPos(int id);
        Eigen::Matrix<double, 3, 4> getFKFeetVel();
        Eigen::Vector3d getFKFeetVel(int i);

        // 一些数学计算函数
        Eigen::Matrix3d quatToRot(const Vector4d& _quat);
        Eigen::Vector3d rotMatToRPY(const Matrix3d& R);
        Eigen::Vector<double, 6> parallelAxis(const Vector3d& dstAxis, const Vector3d& oriAxis, const Vector<double, 6>& oriIm, const double oriM);// 平行轴定理
        Eigen::Matrix3d aI2mI(const Vector<double, 6>& oriIm);
    };

    Body::Body(EstBase* _est,Leg* _legObj[4], double _dt)
    {
        est = _est;
        for (int i = 0; i < 4; i++)
        {
            this->legs[i] = _legObj[i];
        }
        this->dt = _dt;

        Matrix3d Ident = Matrix3d::Identity(); // 初始化一个3维对角矩阵
        dynamicLeft.block<3, 3>(0, 0) = Ident;
        dynamicLeft.block<3, 3>(0, 3) = Ident;
        dynamicLeft.block<3, 3>(0, 6) = Ident;
        dynamicLeft.block<3, 3>(0, 9) = Ident;
    }

    void Body::initParams(Vector3d _leg2body, Vector3d _initRbLegXYPosition, double _Mb[3], Vector<double, 6> _Ib[3], Vector3d _Pb[3])
    {
        this->leg2body[LF] = _leg2body;
        Vector3d tmp = _leg2body;
        tmp(1) = tmp(1) * -1;
        this->leg2body[RF] = tmp;
        tmp(0) = tmp(0) * -1;
        this->leg2body[RB] = tmp;
        tmp(1) = tmp(1) * -1;
        this->leg2body[LB] = tmp;
        this->initRbLegXYPosition = _initRbLegXYPosition;
      /*  this->M = Eigen::Map<Matrix3d>(_M.data(),3,3);
        this->Ib = Eigen::Map<Matrix3d>(_Ib.data(),3,3);*/
        for (int i = 0; i < 3; i++)
        {
            this->Mb[i] = _Mb[i];
            this->Ib[i] = _Ib[i];
            this->Pb[i] = _Pb[i];
        }
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

    void Body::calTbs(int8_t direction)
    {
        if (direction == 1)
        {
            // 当前
            // 三轴欧拉角旋转矩阵
            //Eigen::AngleAxisd rotationx_c(currentBodyState.Ang_xyz(0), Eigen::Vector3d::UnitX());
            //Eigen::AngleAxisd rotationy_c(currentBodyState.Ang_xyz(1), Eigen::Vector3d::UnitY());
            //Eigen::AngleAxisd rotationz_c(currentBodyState.Ang_xyz(2), Eigen::Vector3d::UnitZ());
            //Eigen::Matrix3d rotation_c = (rotationx_c * rotationy_c * rotationz_c).toRotationMatrix();
            //// 构建齐次变换矩阵
            //Rsb_c = rotation_c;
            //Tsb_c.block<3, 3>(0, 0) = rotation_c;
            //Tsb_c.block<3, 1>(0, 3) = currentWorldState.dist;
            // 使用四元数得到变换矩阵
            Eigen::Matrix3d rotation_c = quatToRot(currentBodyState.Quat);
            Rsb_c = rotation_c;
            currentBodyState.Ang_xyz = rotMatToRPY(rotation_c);
            Tsb_c.block<3, 3>(0, 0) = rotation_c;
            Tsb_c.block<3, 1>(0, 3) = currentWorldState.dist;
            // 更新其他世界坐标系下的变量
            currentWorldState.angVel_xyz = rotation_c * currentBodyState.angVel_xyz;
            //currentWorldState.linVel_xyz = rotation_c * currentBodyState.linVel_xyz;
            // currentWorldState.angAcc_xyz = rotation_c * currentBodyState.angAcc_xyz;
            currentWorldState.linAcc_xyz = rotation_c * currentBodyState.linAcc_xyz;
            est->updateTsb(Tsb_c);
        }
        else if (direction == -1)
        {
            // 目标
            // 三轴欧拉角旋转矩阵
            Eigen::AngleAxisd rotationx_t(targetBodyState.Ang_xyz(0), Eigen::Vector3d::UnitX());
            Eigen::AngleAxisd rotationy_t(targetBodyState.Ang_xyz(1), Eigen::Vector3d::UnitY());
            Eigen::AngleAxisd rotationz_t(targetBodyState.Ang_xyz(2), Eigen::Vector3d::UnitZ());
            Eigen::Matrix3d rotation_t = (rotationx_t * rotationy_t * rotationz_t).toRotationMatrix();
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
        for (int i = 0; i < 4; i++)
        {
            // 将机身旋转角速度向量转换成反对称矩阵
            Matrix3d w = v3_to_m3(currentBodyState.angVel_xyz);
            // 更新足端相对于机身坐标系的速度
            currentBodyState.leg_b[i].Velocity = legs[i]->currentLeg.Velocity;
            // 计算足端相对于世界坐标系的速度
            currentWorldState.leg_s[i].Velocity = Rsb_c * (w * currentBodyState.leg_b[i].Position + legs[i]->currentLeg.Velocity);
        }
    }

    void Body::updateDynamic()
    {
        // 首先计算整机重量惯量参数
        double _M = 0;
        Vector3d _P = Vector3d::Zero();
        Vector<double, 6> _I = Vector<double, 6>::Zero();
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                _M += this->legs[i]->Mc[j];
                _P += this->legs[i]->Mc[j] * (this->legs[i]->Pcl[j] + this->leg2body[i]);
            }
        }
        for (int i = 0; i < 3; i++)
        {
            _M += this->Mb[i];
            _P += this->Mb[i] * this->Pb[i];
        }
        _P = _P / _M;// 计算整机重心位置
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                _I += parallelAxis(_P, this->legs[i]->Pcl[j] + this->leg2body[i], this->legs[i]->Ic[j], this->legs[i]->Mc[j]);
            }
        }
        for (int i = 0; i < 3; i++)
        {
            _I += parallelAxis(_P, Pb[i], Ib[i], Mb[i]);
        }
        // 得到整机刚体参数
        this->P = _P;
        this->M = _M * Matrix3d::Identity();
        this->I = aI2mI(_I);

        for (int i = 0; i < 4; i++)
        {
            //dynamicLeft.block<3, 3>(0, i * 3) = Rsb_c;
            Vector3d Pbi = Vector3d::Zero();
            Pbi = Rsb_c * (currentBodyState.leg_b[i].Position - P);
            //Pbi = (currentBodyState.leg_b[i].Position - Pb);
            dynamicLeft.block<3, 3>(3, i * 3) = v3_to_m3(Pbi);
        }
        dynamicRight.block<3, 3>(0, 0) = M;
        dynamicRight.block<3, 3>(3, 3) = Rsb_c * I * Rsb_c.transpose();
        //dynamicRight.block<3, 3>(3, 3) = Ib;
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

    void Body::updateBodyImu(Vector4d _imuQ)
    {
        currentBodyState.Quat = _imuQ;
    }

    void Body::updateBodyGyro(Vector3d _gyro)
    {
        currentBodyState.angVel_xyz = _gyro;
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

    Eigen::Matrix<double, 3, 4> Body::getFKFeetVel()
    {
        Eigen::Matrix<double, 3, 4> out;
        for (int i = 0; i < 4; i++)
        {
            out.block<3, 1>(0, i) = currentWorldState.leg_s[i].Velocity;
        }
        return out;
    }

    Eigen::Vector3d Body::getFKFeetVel(int id)
    {
        Eigen::Vector3d out;
        out = currentWorldState.leg_s[id].Velocity;
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

    Eigen::Vector3d Body::rotMatToRPY(const Matrix3d& R) {
        Vector3d rpy;
        rpy(0) = atan2(R(2, 1), R(2, 2));
        rpy(1) = asin(-R(2, 0));
        rpy(2) = atan2(R(1, 0), R(0, 0));
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

    Eigen::Matrix3d Body::aI2mI(const Vector<double, 6>& oriIm)
    {
        Matrix3d mI;
        mI(0, 0) = oriIm(0);
        mI(1, 1) = oriIm(1);
        mI(2, 2) = oriIm(2);
        mI(0, 1) = oriIm(3);
        mI(1, 0) = oriIm(3);
        mI(0, 2) = oriIm(4);
        mI(2, 0) = oriIm(4);
        mI(1, 2) = oriIm(5);
        mI(2, 1) = oriIm(5);
        return mI;
    }
}
