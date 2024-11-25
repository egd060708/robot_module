#pragma once
#include "Body.h"
#include "mpcCal.h"
#include "PIDmethod.h"
#include "LegCtrl.h"
#include "multiCircle.h"
#include <iostream>

using namespace Eigen;

namespace Quadruped
{
    // 用于平衡控制器的状态
    typedef struct _balanceState
    {
        // 机体状态
        Vector3d p = Vector3d::Zero();
        Vector3d p_dot = Vector3d::Zero();
        Vector3d r = Vector3d::Zero();
        Vector3d r_dot = Vector3d::Zero();
        // 末端执行器状态
        Vector4d pe = Vector4d::Zero();
        Vector4d pe_dot = Vector4d::Zero();
        
    } balanceState;

    class CtrlBase
    {
    protected:

    public:
        // 用于平衡控制器的状态变量
        balanceState targetBalanceState;
        balanceState currentBalanceState;
        // 被控制的机器人本体
        Body* bodyObject;
        // 腿部控制器
        LegCtrl* legsCtrl[4];
        // 动力学模型矩阵
        MatrixXd dynamicLeft;
        MatrixXd dynamicRight;
        // mpc控制器输出
        MatrixXd mpcOut;
        // mpc控制器矩阵
        MatrixXd A;
        MatrixXd B;
        // 控制器权重参数
        MatrixXd Q; // 状态权重矩阵
        MatrixXd F; // 终端补偿矩阵
        MatrixXd R; // 输入权重矩阵
        MatrixXd W; // 输入平滑矩阵
        // 约束矩阵
        MatrixXd lb; // 输入约束
        MatrixXd ub;
        MatrixXd cA; // box约束矩阵
        MatrixXd Alb; // box约束边缘
        MatrixXd Aub;
        // 状态量
        VectorXd y;
        VectorXd x;
        // 其他参数
        double dt;
        bool is_init = false;

        // 构造函数
        CtrlBase(Body* _obj, LegCtrl* _legsCtrl[4], int timeStep)
        {
            bodyObject = _obj;
            for (int i = 0; i < 4; i++)
            {
                this->legsCtrl[i] = _legsCtrl[i];
            }
            this->dt = static_cast<double>(timeStep) * 0.001;
        }
        // 导入权重参数
        virtual void importWeight(const VectorXd& _Q, const VectorXd& _F, const VectorXd& _R, const VectorXd& _W) = 0;
        // 更新动力学模型
        virtual void updateDynamic() = 0;
        // 更新当前状态
        virtual void updateBalanceState() = 0;
        // 执行mpc控制器
        virtual void mpc_adjust(const VectorX<bool>& _enList) = 0;
    };

    class QpCtrl : public CtrlBase {
    private:

    public:
        // 机器人平衡控制器
        mpcCal<6, 12, 20, 1, 5> balanceController;
        double u = 2;// 摩擦系数
        double force_c = 200;// 输出限制

        // 用于机身位控的pid
        PIDmethod linPID[3];
        PIDmethod angPID[3];

        // 构造函数
        QpCtrl(Body* _obj, LegCtrl* _legsCtrl[4], int timeStep):CtrlBase(_obj, _legsCtrl, timeStep)
        {
            dynamicLeft.resize(6, 12);
            dynamicRight.resize(6, 6);
            mpcOut.resize(3, 4);
            A.resize(6, 6);
            B.resize(6, 12);
            Q.resize(6, 6);
            F.resize(6, 6);
            R.resize(12, 12);
            W.resize(12, 12);
            lb.resize(12, 1);
            ub.resize(12, 1);
            cA.resize(20,12);
            Alb.resize(20, 1);
            Aub.resize(20, 1);
            y.resize(6);
            x.resize(6);

            dynamicRight.setZero();
            dynamicLeft.block<3, 3>(0, 0).setIdentity();
            dynamicLeft.block<3, 3>(0, 3).setIdentity();
            dynamicLeft.block<3, 3>(0, 6).setIdentity();
            dynamicLeft.block<3, 3>(0, 9).setIdentity();
            mpcOut.setZero();
            A.setZero();
            B.setZero();
            Q.setZero();
            R.setZero();
            lb.setConstant(-force_c);
            ub.setConstant(force_c);
            cA.setZero();

            // unitree constrain
            Eigen::Matrix<double, 5, 3> _cA;
            _cA.setZero();
            _cA << 1, 0, u, -1, 0, u, 0, 1, u, 0, -1, u, 0, 0, 1;
            cA.block<5, 3>(0, 0) = _cA;
            cA.block<5, 3>(5, 3) = _cA;
            cA.block<5, 3>(10, 6) = _cA;
            cA.block<5, 3>(15, 9) = _cA;
            Alb.setZero();
            Aub.setConstant(100000.);
        }
        // 更新机器人动力学方程（描述为 left*[f] = right）
        void updateDynamic() override;
        // 导入权重参数
        void importWeight(const VectorXd& _Q, const VectorXd& _F, const VectorXd& _R, const VectorXd& _W) override;
        void importPDparam(const Vector<double, 9>& lin,const Vector<double, 9>& ang);
        // 更新当前状态
        void updateBalanceState() override;
        // 设置终端目标
        void setPositionTarget(const Vector3d& _p,const Vector3d& _r);
        void setVelocityTarget(const Vector3d& _p_dot,const Vector3d& _r_dot);
        // 执行mpc控制器
        void mpc_adjust(const VectorX<bool>& _enList) override;
        // 设置腿部接触约束
        void setContactConstrain(const Vector4i& _contact);
        // 直接设置四足输出
        void setLegsForce(const Eigen::Matrix<double, 3, 4>& _force);
    };

    void QpCtrl::updateDynamic()
    {
        for (int i = 0; i < 4; i++)
        {
            //dynamicLeft.block<3, 3>(0, i * 3) = Rsb_c;
            Vector3d Pbi = Vector3d::Zero();
            Pbi = bodyObject->Rsb_c * (bodyObject->currentBodyState.leg_b[i].Position - bodyObject->P);
            //Pbi = (currentBodyState.leg_b[i].Position - Pb);
            dynamicLeft.block<3, 3>(3, i * 3) = bodyObject->v3_to_m3(Pbi);
        }
        dynamicRight.block<3, 3>(0, 0) = bodyObject->M;
        dynamicRight.block<3, 3>(3, 3) = bodyObject->Rsb_c * bodyObject->I * bodyObject->Rsb_c.transpose();
    }

    void QpCtrl::importWeight(const VectorXd& _Q, const VectorXd& _F, const VectorXd& _R, const VectorXd& _W)
    {
        for (int i = 0; i < 6; i++)
        {
            Q(i, i) = _Q(i);
            F(i, i) = _F(i);
        }
        for (int j = 0; j < 12; j++)
        {
            R(j, j) = _R(j);
            W(j, j) = _W(j);
        }
    }

    void QpCtrl::importPDparam(const Vector<double, 9>& lin, const Vector<double, 9>& ang)
    {
        for (int i = 0; i < 3; i++)
        {
            linPID[i].Params_Config(lin(0 + i * 3), 0, lin(1 + i * 3), 0, abs(lin(2 + i * 3)), -abs(lin(2 + i * 3)));
            angPID[i].Params_Config(ang(0 + i * 3), 0, ang(1 + i * 3), 0, abs(ang(2 + i * 3)), -abs(ang(2 + i * 3)));
        }
    }

    void QpCtrl::updateBalanceState()
    {
        currentBalanceState.p = bodyObject->currentWorldState.dist;
        currentBalanceState.p_dot = bodyObject->currentWorldState.linVel_xyz;
        currentBalanceState.r = bodyObject->currentBodyState.Ang_xyz;
        currentBalanceState.r_dot = bodyObject->currentWorldState.angVel_xyz;
    }

    void QpCtrl::setPositionTarget(const Vector3d& _p, const Vector3d& _r)
    {
        targetBalanceState.p = _p;
        targetBalanceState.r = _r;
    }

    void QpCtrl::setVelocityTarget(const Vector3d& _p_dot, const Vector3d& _r_dot)
    {
        targetBalanceState.p_dot = _p_dot;
        targetBalanceState.r_dot = _r_dot;
    }

    void QpCtrl::setContactConstrain(const Vector4i& _contact)
    {
        Eigen::Matrix<double, 5, 3> _cA;
        _cA.setZero();
        Eigen::Matrix<double, 5, 1> _Aub;
        _Aub.setZero();
        // 若该腿不触地，则清零约束矩阵
        for (int i = 0; i < 4; i++)
        {
            if (_contact(i) == 1)
            {
                _cA << 1, 0, u, -1, 0, u, 0, 1, u, 0, -1, u, 0, 0, 1;
                _Aub.setConstant(100000.);
            }
            else
            {
                _cA << 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                _Aub.setZero();
            }
            this->cA.block<5, 3>(5 * i, 3 * i) = _cA;
            this->Aub.block<5, 1>(5 * i, 0) = _Aub;
        }
    }

    void QpCtrl::mpc_adjust(const VectorX<bool>& _enList)
    {
        static multiCircle angC[3] = { multiCircle(3.1415926),multiCircle(3.1415926),multiCircle(3.1415926) };
        for (int i = 0; i < 3; i++)
        {
            if (_enList(2 * i) == true)
            {
                linPID[i].target = targetBalanceState.p(i);
                linPID[i].current = currentBalanceState.p(i);
                linPID[i].Adjust(0, currentBalanceState.p_dot(i));
                targetBalanceState.p_dot(i) = linPID[i].out;
            }
            else
            {
                //targetBalanceState.p_dot(i) = currentBalanceState.p_dot(i);
            }

            if (_enList(2 * i + 1) == true)
            {
                angPID[i].target = targetBalanceState.r(i);
                angPID[i].current = angC[i].f(currentBalanceState.r(i));
                angPID[i].Adjust(0, bodyObject->currentBodyState.angVel_xyz(i));
                targetBalanceState.r_dot(i) = angPID[i].out;
                if (i == 2)
                {
                    targetBalanceState.r_dot = bodyObject->Rsb_c * targetBalanceState.r_dot;
                }
            }
            else
            {
                //targetBalanceState.r_dot(i) = currentBalanceState.r_dot(i);
            }
        }

        B = dynamicRight.inverse() * dynamicLeft;
        balanceController.setConstrain(lb, ub);
        balanceController.setBoxConstrain(cA, Alb, Aub);
        y.block<3, 1>(0, 0) = targetBalanceState.p_dot;
        y.block<3, 1>(3, 0) = targetBalanceState.r_dot;
        x.block<3, 1>(0, 0) = currentBalanceState.p_dot;
        x.block<3, 1>(3, 0) = currentBalanceState.r_dot;
        balanceController.mpc_update(y, x, 100, 0.002);
        balanceController.mpc_init(A, B, Q, F, R, W, dt);
        balanceController.mpc_solve();
        for (int i = 0; i < 4; i++)
        {
            this->mpcOut.col(i) = -bodyObject->Rsb_c.transpose() * balanceController.getOutput().block<3, 1>(3 * i, 0);
            //this->mpcOut.col(i) = -balanceController.getOutput().block<3, 1>(3 * i, 0);
        }
    }

    void QpCtrl::setLegsForce(const Eigen::Matrix<double, 3, 4>& _force)
    {
        for (int i = 0; i < 4; i++)
        {
            this->legsCtrl[i]->setEndForceTar(_force.col(i));
        }
    }

    class QpwCtrl : public CtrlBase {
    private:

    public:
        // 机器人平衡控制器
        mpcCal<10, 16, 20, 1, 5> balanceController;
        double u = 2;// 摩擦系数
        double force_c = 200;// 输出力限制限制
        double tau_c = 20; // 输出力矩限制

        // 用于机身位控的pid
        PIDmethod linPID[3];
        PIDmethod angPID[3];

        // 用于轮速控制的pid
        PIDmethod wheelPID[4];

        // 构造函数
        QpwCtrl(Body* _obj, LegCtrl* _legsCtrl[4], int timeStep) :CtrlBase(_obj, _legsCtrl, timeStep)
        {
            dynamicLeft.resize(10, 16);
            dynamicRight.resize(10, 10);
            mpcOut.resize(4, 4);
            A.resize(10, 10);
            B.resize(10, 16);
            Q.resize(10, 10);
            F.resize(10, 10);
            R.resize(16, 16);
            W.resize(16, 16);
            lb.resize(16, 1);
            ub.resize(16, 1);
            cA.resize(20, 16);
            Alb.resize(20, 1);
            Aub.resize(20, 1);
            y.resize(10);
            x.resize(10);

            dynamicRight.setZero();
            dynamicRight.block(6, 6, 4, 4) = (bodyObject->legs[0]->Ic[3](1) / pow(bodyObject->legs[0]->L4 + bodyObject->legs[0]->L4b, 2) + bodyObject->legs[0]->Mc[3]) * Matrix4d::Identity();
            dynamicLeft.setZero();
            dynamicLeft.block<3, 3>(0, 0).setIdentity();
            dynamicLeft.block<3, 3>(0, 3).setIdentity();
            dynamicLeft.block<3, 3>(0, 6).setIdentity();
            dynamicLeft.block<3, 3>(0, 9).setIdentity();
            dynamicLeft(6, 0) = -1;
            dynamicLeft(7, 3) = -1;
            dynamicLeft(8, 6) = -1;
            dynamicLeft(9, 9) = -1;
            dynamicLeft.block(6, 12, 4, 4) = (1. / (bodyObject->legs[0]->L4 + bodyObject->legs[0]->L4b)) * Matrix4d::Identity();
            mpcOut.setZero();
            A.setZero();
            B.setZero();
            Q.setZero();
            R.setZero();
            lb.block(0, 0, 12, 1).setConstant(-force_c);
            ub.block(0, 0, 12, 1).setConstant(force_c);
            lb.block(12, 0, 4, 1).setConstant(-tau_c);
            ub.block(12, 0, 4, 1).setConstant(tau_c);
            cA.setZero();

            // unitree constrain
            Eigen::Matrix<double, 5, 3> _cA;
            _cA.setZero();
            _cA << 1, 0, u, -1, 0, u, 0, 1, u, 0, -1, u, 0, 0, 1;
            cA.block<5, 3>(0, 0) = _cA;
            cA.block<5, 3>(5, 3) = _cA;
            cA.block<5, 3>(10, 6) = _cA;
            cA.block<5, 3>(15, 9) = _cA;
            Alb.setZero();
            Aub.setConstant(100000.);
        }

        // 更新机器人动力学方程（描述为 left*[f] = right）
        void updateDynamic() override;
        // 导入权重参数
        void importWeight(const VectorXd& _Q, const VectorXd& _F, const VectorXd& _R, const VectorXd& _W) override;
        void importPDparam(const Vector<double, 9>& lin, const Vector<double, 9>& ang, const Vector<double, 3>& wheel);
        // 更新当前状态
        void updateBalanceState() override;
        // 设置终端目标
        void setPositionTarget(const Vector3d& _p, const Vector3d& _r, const Vector4d& _pe);
        void setVelocityTarget(const Vector3d& _p_dot, const Vector3d& _r_dot, const Vector4d& _pe_dot);
        // 执行mpc控制器
        void mpc_adjust(const VectorX<bool>& _enList) override;
        // 设置接触约束
        void setContactConstrain(const Vector4i& _contact);
        // 直接设置四足输出
        void setLegsForce(const Eigen::Matrix<double, 3, 4>& _force, const Eigen::Vector4d& _tau);
    };

    void QpwCtrl::updateDynamic()
    {
        for (int i = 0; i < 4; i++)
        {
            //dynamicLeft.block<3, 3>(0, i * 3) = Rsb_c;
            Vector3d Pbi = Vector3d::Zero();
            Pbi = bodyObject->Rsb_c * (bodyObject->currentBodyState.leg_b[i].Position - bodyObject->P);
            //Pbi = (currentBodyState.leg_b[i].Position - Pb);
            dynamicLeft.block<3, 3>(3, i * 3) = bodyObject->v3_to_m3(Pbi);
        }
        dynamicRight.block<3, 3>(0, 0) = bodyObject->M;
        dynamicRight.block<3, 3>(3, 3) = bodyObject->Rsb_c * bodyObject->I * bodyObject->Rsb_c.transpose();
        Eigen::Matrix<double, 4, 3> s = Eigen::Matrix<double, 4, 3>::Zero();
        s(0,0) = -1;
        dynamicLeft.block(6, 0, 4, 3) = s * bodyObject->Rsb_c.transpose();
        s(0, 0) = 0;
        s(1, 0) = -1;
        dynamicLeft.block(6, 3, 4, 3) = s * bodyObject->Rsb_c.transpose();
        s(1, 0) = 0;
        s(2, 0) = -1;
        dynamicLeft.block(6, 6, 4, 3) = s * bodyObject->Rsb_c.transpose();
        s(2, 0) = 0;
        s(3, 0) = -1;
        dynamicLeft.block(6, 9, 4, 3) = s * bodyObject->Rsb_c.transpose();
    }

    void QpwCtrl::importWeight(const VectorXd& _Q, const VectorXd& _F, const VectorXd& _R, const VectorXd& _W)
    {
        for (int i = 0; i < 10; i++)
        {
            Q(i, i) = _Q(i);
            F(i, i) = _F(i);
        }
        for (int j = 0; j < 16; j++)
        {
            R(j, j) = _R(j);
            W(j, j) = _W(j);
        }
    }

    void QpwCtrl::importPDparam(const Vector<double, 9>& lin, const Vector<double, 9>& ang, const Vector<double, 3>& wheel)
    {
        for (int i = 0; i < 3; i++)
        {
            linPID[i].Params_Config(lin(0 + i * 3), 0, lin(1 + i * 3), 0, abs(lin(2 + i * 3)), -abs(lin(2 + i * 3)));
            angPID[i].Params_Config(ang(0 + i * 3), 0, ang(1 + i * 3), 0, abs(ang(2 + i * 3)), -abs(ang(2 + i * 3)));
        }
        for (int i = 0; i < 4; i++)
        {
            wheelPID[i].Params_Config(wheel(0), wheel(1), abs(wheel(2)), -abs(wheel(2)));
        }
    }

    void QpwCtrl::updateBalanceState()
    {
        currentBalanceState.p = bodyObject->currentWorldState.dist;
        currentBalanceState.p_dot = bodyObject->currentWorldState.linVel_xyz;
        currentBalanceState.r = bodyObject->currentBodyState.Ang_xyz;
        currentBalanceState.r_dot = bodyObject->currentWorldState.angVel_xyz;
        for (int i = 0; i < 4; i++)
        {
            /*currentBalanceState.pe(i) = bodyObject->est->getEstFeetPosS(i)(0);*/
            currentBalanceState.pe(i) = bodyObject->currentBodyState.leg_b[i].Position(0);
            currentBalanceState.pe_dot(i) = bodyObject->currentBodyState.leg_b[i].VelocityW(0);
        }
        //std::cout << "cbpe: \n" << currentBalanceState.pe_dot << std::endl;
    }

    void QpwCtrl::setPositionTarget(const Vector3d& _p, const Vector3d& _r, const Vector4d& _pe)
    {
        targetBalanceState.p = _p;
        targetBalanceState.r = _r;
        targetBalanceState.pe = _pe;
    }

    void QpwCtrl::setVelocityTarget(const Vector3d& _p_dot, const Vector3d& _r_dot, const Vector4d& _pe_dot)
    {
        targetBalanceState.p_dot = _p_dot;
        targetBalanceState.r_dot = _r_dot;
        targetBalanceState.pe_dot = _pe_dot;
    }

    void QpwCtrl::mpc_adjust(const VectorX<bool>& _enList)
    {
        static multiCircle angC[3] = { multiCircle(3.1415926),multiCircle(3.1415926),multiCircle(3.1415926) };
        /*std::cout << "dLeft :" << dynamicLeft << std::endl;
        std::cout << "dRight:" << dynamicRight << std::endl;
        std::cout << "lb:" << lb << std::endl;
        std::cout << "ub:" << ub << std::endl;
        std::cout << "cA:" << cA << std::endl;
        std::cout << "Alb:" << Alb << std::endl;
        std::cout << "Aub:" << Aub << std::endl;*/
        for (int i = 0; i < 3; i++)
        {
            if (_enList(2 * i) == true)
            {
                linPID[i].target = targetBalanceState.p(i);
                linPID[i].current = currentBalanceState.p(i);
                linPID[i].Adjust(0, currentBalanceState.p_dot(i));
                targetBalanceState.p_dot(i) = linPID[i].out;
            }
            else
            {
                targetBalanceState.p_dot(i) = currentBalanceState.p_dot(i);
            }

            if (_enList(2 * i + 1) == true)
            {
                angPID[i].target = targetBalanceState.r(i);
                angPID[i].current = angC[i].f(currentBalanceState.r(i));
                angPID[i].Adjust(0, bodyObject->currentBodyState.angVel_xyz(i));
                targetBalanceState.r_dot(i) = angPID[i].out;
                if (i == 2)
                {
                    targetBalanceState.r_dot = bodyObject->Rsb_c * targetBalanceState.r_dot;
                }
            }
            else
            {
                targetBalanceState.r_dot(i) = currentBalanceState.r_dot(i);
            }
        }
        for (int i = 0; i < 4; i++)
        {
            wheelPID[i].target = targetBalanceState.pe(i);
            wheelPID[i].current = currentBalanceState.pe(i);
            //std::cout << "wpid: " << wheelPID[i].target << ", " << wheelPID[i].current << std::endl;
            wheelPID[i].Adjust(0);
            targetBalanceState.pe_dot(i) = wheelPID[i].out;
        }

        B = dynamicRight.inverse() * dynamicLeft;
        balanceController.setConstrain(lb, ub);
        balanceController.setBoxConstrain(cA, Alb, Aub);
        y.block<3, 1>(0, 0) = targetBalanceState.p_dot;
        y.block<3, 1>(3, 0) = targetBalanceState.r_dot;
        y.block<4, 1>(6, 0) = targetBalanceState.pe_dot;
        x.block<3, 1>(0, 0) = currentBalanceState.p_dot;
        x.block<3, 1>(3, 0) = currentBalanceState.r_dot;
        x.block<4, 1>(6, 0) = currentBalanceState.pe_dot;
        balanceController.mpc_update(y, x, 100, 0.002);
        balanceController.mpc_init(A, B, Q, F, R, W, dt);
        balanceController.mpc_solve();
        for (int i = 0; i < 4; i++)
        {
            this->mpcOut.block(0,i,3,1) = -bodyObject->Rsb_c.transpose() * balanceController.getOutput().block<3, 1>(3 * i, 0);
            this->mpcOut(3, i) = balanceController.getOutput()(12 + i, 0);
            //std::cout << "wheelPID: \n" << wheelPID[i].target << wheelPID[i].current << wheelPID[i].out << std::endl;
        }
    }

    void QpwCtrl::setContactConstrain(const Vector4i& _contact)
    {
        Eigen::Matrix<double, 5, 3> _cA;
        _cA.setZero();
        Eigen::Matrix<double, 5, 1> _Aub;
        _Aub.setZero();
        // 若该腿不触地，则清零约束矩阵
        for (int i = 0; i < 4; i++)
        {
            if (_contact(i) == 1)
            {
                _cA << 1, 0, u, -1, 0, u, 0, 1, u, 0, -1, u, 0, 0, 1;
                _Aub.setConstant(100000.);
            }
            else
            {
                _cA << 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                _Aub.setZero();
            }
            this->cA.block<5, 3>(5 * i, 3 * i) = _cA;
            this->Aub.block<5, 1>(5 * i, 0) = _Aub;
        }
    }

    void QpwCtrl::setLegsForce(const Eigen::Matrix<double, 3, 4>& _force, const Eigen::Vector4d& _tau)
    {
        for (int i = 0; i < 4; i++)
        {
            this->legsCtrl[i]->setEndForceTar(_force.col(i));
            this->legsCtrl[i]->setEndTauTar(_tau(i));
        }
    }
}
