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
            this->dt = static_cast<double>(timeStep) * 0.005;
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

    class QpPIDVCtrl : public CtrlBase {
    private:

    public:
        // 机器人平衡控制器
        mpcCal<6, 12, 20, 1, 10> balanceController;
        double u = 2;// 摩擦系数
        double force_c = 200;// 输出限制

        // 用于机身位控的pid
        PIDmethod linPID[3];
        PIDmethod angPID[3];

        // 构造函数
        QpPIDVCtrl(Body* _obj, LegCtrl* _legsCtrl[4], int timeStep):CtrlBase(_obj, _legsCtrl, timeStep), balanceController(PL_NONE)
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

    void QpPIDVCtrl::updateDynamic()
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

    void QpPIDVCtrl::importWeight(const VectorXd& _Q, const VectorXd& _F, const VectorXd& _R, const VectorXd& _W)
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

    void QpPIDVCtrl::importPDparam(const Vector<double, 9>& lin, const Vector<double, 9>& ang)
    {
        for (int i = 0; i < 3; i++)
        {
            linPID[i].Params_Config(lin(0 + i * 3), 0, lin(1 + i * 3), 0, abs(lin(2 + i * 3)), -abs(lin(2 + i * 3)));
            angPID[i].Params_Config(ang(0 + i * 3), 0, ang(1 + i * 3), 0, abs(ang(2 + i * 3)), -abs(ang(2 + i * 3)));
        }
    }

    void QpPIDVCtrl::updateBalanceState()
    {
        currentBalanceState.p = bodyObject->currentWorldState.dist;
        currentBalanceState.p_dot = bodyObject->currentWorldState.linVel_xyz;
        currentBalanceState.r = bodyObject->currentBodyState.Ang_xyz;
        currentBalanceState.r_dot = bodyObject->currentWorldState.angVel_xyz;
    }

    void QpPIDVCtrl::setPositionTarget(const Vector3d& _p, const Vector3d& _r)
    {
        targetBalanceState.p = _p;
        targetBalanceState.r = _r;
    }

    void QpPIDVCtrl::setVelocityTarget(const Vector3d& _p_dot, const Vector3d& _r_dot)
    {
        targetBalanceState.p_dot = _p_dot;
        targetBalanceState.r_dot = _r_dot;
    }

    void QpPIDVCtrl::setContactConstrain(const Vector4i& _contact)
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
        lb.setConstant(-force_c);
        ub.setConstant(force_c);
    }

    void QpPIDVCtrl::mpc_adjust(const VectorX<bool>& _enList)
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
        balanceController.mpc_solve(0);
        for (int i = 0; i < 4; i++)
        {
            this->mpcOut.col(i) = -bodyObject->Rsb_c.transpose() * balanceController.getOutput().block<3, 1>(3 * i, 0);
            //this->mpcOut.col(i) = -balanceController.getOutput().block<3, 1>(3 * i, 0);
        }
    }

    void QpPIDVCtrl::setLegsForce(const Eigen::Matrix<double, 3, 4>& _force)
    {
        for (int i = 0; i < 4; i++)
        {
            this->legsCtrl[i]->setEndForceTar(_force.col(i));
        }
    }


    class QpPVCtrl : public CtrlBase {
    private:

    public:
        // 机器人平衡控制器
        mpcCal<15, 12, 20, 1, 5> balanceController;
        double u = 0.8;// 摩擦系数
        double force_c = 200;// 输出限制
        Eigen::Vector3d g = Eigen::Vector3d(0, 0, -9.81);
        PIDmethod wheelPID[4];
        Eigen::Vector4d wheelTau;

        // 构造函数
        QpPVCtrl(Body* _obj, LegCtrl* _legsCtrl[4], int timeStep) :CtrlBase(_obj, _legsCtrl, timeStep), balanceController(PL_NONE)
        {
            dynamicLeft.resize(6, 12);
            dynamicRight.resize(6, 6);
            mpcOut.resize(3, 4);
            A.resize(15, 15);
            B.resize(15, 12);
            Q.resize(15, 15);
            F.resize(15, 15);
            R.resize(12, 12);
            W.resize(12, 12);
            lb.resize(12, 1);
            ub.resize(12, 1);
            cA.resize(20, 12);
            Alb.resize(20, 1);
            Aub.resize(20, 1);
            y.resize(15);
            x.resize(15);

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
        // 更新当前状态
        void updateBalanceState() override;
        // 设置终端目标
        void setPositionTarget(const Vector3d& _p, const Vector3d& _r);
        void setVelocityTarget(const Vector3d& _p_dot, const Vector3d& _r_dot);
        // 执行mpc控制器
        void mpc_adjust(const VectorX<bool>& _enList) override;
        // 执行轮速控制器
        void wheel_adjust(const Vector4d& _p, const Vector4d& _p_dot);
        // 设置腿部接触约束
        void setContactConstrain(const Vector4i& _contact, const Matrix3d& _extR);
        // 直接设置四足输出
        void setLegsForce(const Eigen::Matrix<double, 3, 4>& _force);
        void setLegsForce(const Eigen::Matrix<double, 3, 4>& _force, const Eigen::Vector4d& _tau);
    };

    void QpPVCtrl::updateDynamic()
    {
        /*Eigen::Vector3d simpleP = 0.1 * bodyObject->P;
        Eigen::Matrix3d simpleI = Eigen::Vector3d(bodyObject->I(0, 0), bodyObject->I(1, 1), bodyObject->I(2, 2)).asDiagonal();*/
        for (int i = 0; i < 4; i++)
        {
            Vector3d Pbi = Vector3d::Zero();
            Pbi = bodyObject->Rsb_c * (bodyObject->currentBodyState.leg_b[i].Position - bodyObject->P);
            //Pbi = bodyObject->Rsb_c * (bodyObject->currentBodyState.leg_b[i].Position - simpleP);
            dynamicLeft.block<3, 3>(3, i * 3) = bodyObject->v3_to_m3(Pbi);
        }
        dynamicRight.block<3, 3>(0, 0) = bodyObject->M;
        dynamicRight.block<3, 3>(3, 3) = bodyObject->Rsb_c * bodyObject->I * bodyObject->Rsb_c.transpose();
        //dynamicRight.block<3, 3>(3, 3) = bodyObject->Rsb_c * simpleI * bodyObject->Rsb_c.transpose();
        /*std::cout << "M:" << bodyObject->M << std::endl;
        std::cout << "I:" << simpleI << std::endl;
        std::cout << "P:" << simpleP << std::endl;*/
    }

    void QpPVCtrl::importWeight(const VectorXd& _Q, const VectorXd& _F, const VectorXd& _R, const VectorXd& _W)
    {
        for (int i = 0; i < 15; i++)
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

    void QpPVCtrl::updateBalanceState()
    {
        static multiCircle angC = multiCircle(3.1415926);
        currentBalanceState.p = bodyObject->currentWorldState.dist;
        currentBalanceState.p_dot = bodyObject->currentWorldState.linVel_xyz;
        currentBalanceState.r(0) = bodyObject->currentBodyState.Ang_xyz(0);
        currentBalanceState.r(1) = bodyObject->currentBodyState.Ang_xyz(1);
        currentBalanceState.r(2) = angC.f(bodyObject->currentBodyState.Ang_xyz(2));
        currentBalanceState.r_dot = bodyObject->currentWorldState.angVel_xyz;
    }

    void QpPVCtrl::setPositionTarget(const Vector3d& _p, const Vector3d& _r)
    {
        targetBalanceState.p = _p;
        targetBalanceState.r = _r;
    }

    void QpPVCtrl::setVelocityTarget(const Vector3d& _p_dot, const Vector3d& _r_dot)
    {
        targetBalanceState.p_dot = _p_dot;
        targetBalanceState.r_dot = _r_dot;
    }

    void QpPVCtrl::setContactConstrain(const Vector4i& _contact, const Matrix3d& _extR)
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
                /*_cA << _extR(0, 0) + u * _extR(2, 0), _extR(0, 1) + u * _extR(2, 1), _extR(0, 2) + u * _extR(2, 2), \
                    - _extR(0, 0) + u * _extR(2, 0), -_extR(0, 1) + u * _extR(2, 1), -_extR(0, 2) + u * _extR(2, 2), \
                    _extR(1, 0) + u * _extR(2, 0), _extR(1, 1) + u * _extR(2, 1), _extR(1, 2) + u * _extR(2, 2), \
                    - _extR(1, 0) + u * _extR(2, 0), -_extR(1, 1) + u * _extR(2, 1), -_extR(1, 2) + u * _extR(2, 2), \
                    _extR(2, 0), _extR(2, 1), _extR(2, 2);*/
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
        lb.setConstant(-force_c);
        ub.setConstant(force_c);
    }

    void QpPVCtrl::mpc_adjust(const VectorX<bool>& _enList)
    {
        A.block(0, 6, 3, 3) = Eigen::Matrix3d::Identity();
        //A.block(3, 9, 3, 3) = Eigen::Matrix3d::Identity();
        A.block(3, 9, 3, 3) = bodyObject->Rsb_c.transpose();
        A.block(6, 12, 3, 3) = Eigen::Matrix3d::Identity();
        B.block(6,0,6,12) = dynamicRight.inverse() * dynamicLeft;
        balanceController.setConstrain(lb, ub);
        balanceController.setBoxConstrain(cA, Alb, Aub);
        y.block<3, 1>(0, 0) = targetBalanceState.p;
        y.block<3, 1>(3, 0) = targetBalanceState.r;
        y.block<3, 1>(6, 0) = targetBalanceState.p_dot;
        y.block<3, 1>(9, 0) = targetBalanceState.r_dot;
        y.block<3, 1>(12, 0) = g;
        x.block<3, 1>(0, 0) = currentBalanceState.p;
        x.block<3, 1>(3, 0) = currentBalanceState.r;
        x.block<3, 1>(6, 0) = currentBalanceState.p_dot;
        x.block<3, 1>(9, 0) = currentBalanceState.r_dot;
        x.block<3, 1>(12, 0) = g;
        /*std::cout << "targetP: " << targetBalanceState.p << std::endl;
        std::cout << "currenP: " << currentBalanceState.p << std::endl;*/
        /*std::cout << "targetPdot: " << targetBalanceState.p_dot << std::endl;
        std::cout << "currenPdot: " << currentBalanceState.p_dot << std::endl;*/
        /*std::cout << "targetR: " << targetBalanceState.r << std::endl;
        std::cout << "currenR: " << currentBalanceState.r << std::endl;
        std::cout << "targetRdot: " << targetBalanceState.r_dot << std::endl;
        std::cout << "currenRdot: " << currentBalanceState.r_dot << std::endl;*/
        balanceController.mpc_update(y, x, 100, 0.02);
        balanceController.mpc_init(A, B, Q, F, R, W, dt);
        balanceController.mpc_solve(0);
        for (int i = 0; i < 4; i++)
        {
            this->mpcOut.col(i) = -bodyObject->Rsb_c.transpose() * balanceController.getOutput().block<3, 1>(3 * i, 0);
            //this->mpcOut.col(i) = -balanceController.getOutput().block<3, 1>(3 * i, 0);
        }
    }

    void QpPVCtrl::wheel_adjust(const Vector4d& _p, const Vector4d& _p_dot)
    {
        for (int i = 0; i < 4; i++)
        {
            wheelPID[i].Params_Config(200, 0., -0.2, 0.5, 20, -20);
            wheelPID[i].target = _p(i);
            wheelPID[i].current = bodyObject->currentBodyState.leg_b[i].Position(0);
            double feedforward = (bodyObject->Rsb_c.transpose() * balanceController.getOutput().block<3, 1>(3 * i, 0))(0) * bodyObject->legs[i]->Reff;
            wheelTau(i) = feedforward + wheelPID[i].Adjust(0, bodyObject->currentBodyState.leg_b[i].VelocityW(0));
        }
    }

    void QpPVCtrl::setLegsForce(const Eigen::Matrix<double, 3, 4>& _force)
    {
        for (int i = 0; i < 4; i++)
        {
            this->legsCtrl[i]->setEndForceTar(_force.col(i));
        }
    }

    void QpPVCtrl::setLegsForce(const Eigen::Matrix<double, 3, 4>& _force, const Eigen::Vector4d& _tau)
    {
        for (int i = 0; i < 4; i++)
        {
            this->legsCtrl[i]->setEndForceTar(_force.col(i));
            this->legsCtrl[i]->setEndTauTar(_tau(i));
        }
    }

    class QpwPIDVCtrl : public CtrlBase {
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
        QpwPIDVCtrl(Body* _obj, LegCtrl* _legsCtrl[4], int timeStep) :CtrlBase(_obj, _legsCtrl, timeStep), balanceController(PL_NONE)
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
            //dynamicRight.block(6, 6, 4, 4) = (bodyObject->legs[0]->Ic[3](1) / pow(bodyObject->legs[0]->L4 + bodyObject->legs[0]->L4b, 2) + bodyObject->legs[0]->Mc[3]) * Matrix4d::Identity();
            dynamicLeft.setZero();
            dynamicLeft.block<3, 3>(0, 0).setIdentity();
            dynamicLeft.block<3, 3>(0, 3).setIdentity();
            dynamicLeft.block<3, 3>(0, 6).setIdentity();
            dynamicLeft.block<3, 3>(0, 9).setIdentity();
            dynamicLeft(6, 0) = -1;
            dynamicLeft(7, 3) = -1;
            dynamicLeft(8, 6) = -1;
            dynamicLeft(9, 9) = -1;
            //dynamicLeft.block(6, 12, 4, 4) = (1. / (bodyObject->legs[0]->L4 + bodyObject->legs[0]->L4b)) * Matrix4d::Identity();
            for (int i = 0; i < 4; i++)
            {
                dynamicRight(6 + i, 6 + i) = bodyObject->legs[i]->Ic[3](1) / pow(bodyObject->legs[i]->Reff, 2) + bodyObject->legs[i]->Mc[3];
                dynamicLeft(6 + i, 12 + i) = 1. / bodyObject->legs[i]->Reff;
            }
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

    void QpwPIDVCtrl::updateDynamic()
    {
        for (int i = 0; i < 4; i++)
        {
            // 四足接触点位置反对称矩阵的计算
            //dynamicLeft.block<3, 3>(0, i * 3) = Rsb_c;
            Vector3d Pbi = Vector3d::Zero();
            Pbi = bodyObject->Rsb_c * (bodyObject->currentBodyState.leg_b[i].Position - bodyObject->P);
            //Pbi = (currentBodyState.leg_b[i].Position - Pb);
            dynamicLeft.block<3, 3>(3, i * 3) = bodyObject->v3_to_m3(Pbi);

            // 轮相关物理参数更新
            dynamicRight(6 + i, 6 + i) = bodyObject->legs[i]->Ic[3](1) / pow(bodyObject->legs[i]->Reff, 2) + bodyObject->legs[i]->Mc[3];
            dynamicLeft(6 + i, 12 + i) = 1. / bodyObject->legs[i]->Reff;
        }
        dynamicRight.block<3, 3>(0, 0) = bodyObject->M;
        dynamicRight.block<3, 3>(3, 3) = bodyObject->Rsb_c * bodyObject->I * bodyObject->Rsb_c.transpose();
        
        // 为了减少计算量，轮速规划是基于车体坐标系的，因此要对世界坐标系的力做映射
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

    void QpwPIDVCtrl::importWeight(const VectorXd& _Q, const VectorXd& _F, const VectorXd& _R, const VectorXd& _W)
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

    void QpwPIDVCtrl::importPDparam(const Vector<double, 9>& lin, const Vector<double, 9>& ang, const Vector<double, 3>& wheel)
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

    void QpwPIDVCtrl::updateBalanceState()
    {
        currentBalanceState.p = bodyObject->currentWorldState.dist;
        currentBalanceState.p_dot = bodyObject->currentWorldState.linVel_xyz;
        currentBalanceState.r = bodyObject->currentBodyState.Ang_xyz;
        currentBalanceState.r_dot = bodyObject->currentWorldState.angVel_xyz;
        for (int i = 0; i < 4; i++)
        {
            /*currentBalanceState.pe(i) = bodyObject->est->getEstFeetPosB(i)(0);
            currentBalanceState.pe_dot(i) = bodyObject->est->getEstFeetVelB(i)(0);*/
            currentBalanceState.pe(i) = bodyObject->currentBodyState.leg_b[i].Position(0);
            currentBalanceState.pe_dot(i) = bodyObject->currentBodyState.leg_b[i].VelocityW(0);
        }
        //std::cout << "cbpe: \n" << currentBalanceState.pe_dot << std::endl;
    }

    void QpwPIDVCtrl::setPositionTarget(const Vector3d& _p, const Vector3d& _r, const Vector4d& _pe)
    {
        targetBalanceState.p = _p;
        targetBalanceState.r = _r;
        targetBalanceState.pe = _pe;
    }

    void QpwPIDVCtrl::setVelocityTarget(const Vector3d& _p_dot, const Vector3d& _r_dot, const Vector4d& _pe_dot)
    {
        targetBalanceState.p_dot = _p_dot;
        targetBalanceState.r_dot = _r_dot;
        targetBalanceState.pe_dot = _pe_dot;
    }

    void QpwPIDVCtrl::mpc_adjust(const VectorX<bool>& _enList)
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
        balanceController.mpc_solve(0);
        for (int i = 0; i < 4; i++)
        {
            this->mpcOut.block(0,i,3,1) = -bodyObject->Rsb_c.transpose() * balanceController.getOutput().block<3, 1>(3 * i, 0);
            this->mpcOut(3, i) = balanceController.getOutput()(12 + i, 0);
            //std::cout << "wheelPID: \n" << wheelPID[i].target << wheelPID[i].current << wheelPID[i].out << std::endl;
        }
    }

    void QpwPIDVCtrl::setContactConstrain(const Vector4i& _contact)
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
        lb.block(0, 0, 12, 1).setConstant(-force_c);
        ub.block(0, 0, 12, 1).setConstant(force_c);
        lb.block(12, 0, 4, 1).setConstant(-tau_c);
        ub.block(12, 0, 4, 1).setConstant(tau_c);
    }

    void QpwPIDVCtrl::setLegsForce(const Eigen::Matrix<double, 3, 4>& _force, const Eigen::Vector4d& _tau)
    {
        for (int i = 0; i < 4; i++)
        {
            this->legsCtrl[i]->setEndForceTar(_force.col(i));
            this->legsCtrl[i]->setEndTauTar(_tau(i));
        }
    }


    class QpwPVCtrl : public CtrlBase {
    private:

    public:
        // 机器人平衡控制器
        mpcCal<23, 16, 28, 1, 5> balanceController;
        double u = 0.8;// 摩擦系数
        double force_cz = 750;// 输出竖直力限制限制
        double force_cxy = 750;// 输出水平力矩限制
        double tau_c = 20; // 输出力矩限制
        Eigen::Vector3d g = Eigen::Vector3d(0, 0, -9.81);
        double fftauRatio[4] = { 0.5 };

        // 构造函数
        QpwPVCtrl(Body* _obj, LegCtrl* _legsCtrl[4], int timeStep) :CtrlBase(_obj, _legsCtrl, timeStep), balanceController(PL_LOW)
        {
            dynamicLeft.resize(10, 16);
            dynamicRight.resize(10, 10);
            mpcOut.resize(4, 4);
            A.resize(23, 23);
            B.resize(23, 16);
            Q.resize(23, 23);
            F.resize(23, 23);
            R.resize(16, 16);
            W.resize(16, 16);
            lb.resize(16, 1);
            ub.resize(16, 1);
            cA.resize(28, 16);
            Alb.resize(28, 1);
            Aub.resize(28, 1);
            y.resize(23);
            x.resize(23);

            dynamicRight.setZero();
            //dynamicRight.block(6, 6, 4, 4) = (bodyObject->legs[0]->Ic[3](1) / pow(bodyObject->legs[0]->L4 + bodyObject->legs[0]->L4b, 2) + bodyObject->legs[0]->Mc[3]) * Matrix4d::Identity();
            dynamicLeft.setZero();
            dynamicLeft.block<3, 3>(0, 0).setIdentity();
            dynamicLeft.block<3, 3>(0, 3).setIdentity();
            dynamicLeft.block<3, 3>(0, 6).setIdentity();
            dynamicLeft.block<3, 3>(0, 9).setIdentity();
            dynamicLeft(6, 0) = -1;
            dynamicLeft(7, 3) = -1;
            dynamicLeft(8, 6) = -1;
            dynamicLeft(9, 9) = -1;
            //dynamicLeft.block(6, 12, 4, 4) = (1. / (bodyObject->legs[0]->L4 + bodyObject->legs[0]->L4b)) * Matrix4d::Identity();
            for (int i = 0; i < 4; i++)
            {
                dynamicRight(6 + i, 6 + i) = bodyObject->legs[i]->Ic[3](1) / pow(bodyObject->legs[i]->Reff, 2) + bodyObject->legs[i]->Mc[3];
                dynamicLeft(6 + i, 12 + i) = 1. / bodyObject->legs[i]->Reff;
            }
            mpcOut.setZero();
            A.setZero();
            B.setZero();
            Q.setZero();
            R.setZero();
            for (int i = 0; i < 4; i++)
            {
                lb(3 * i, 0) = -force_cxy;
                lb(3 * i + 1, 0) = -force_cxy;
                lb(3 * i + 2, 0) = -force_cz;
                ub(3 * i, 0) = force_cxy;
                ub(3 * i + 1, 0) = force_cxy;
                ub(3 * i + 2, 0) = force_cz;
            }
            lb.block(12, 0, 4, 1).setConstant(-tau_c);
            ub.block(12, 0, 4, 1).setConstant(tau_c);
            cA.setZero();

            // unitree constrain
            Eigen::Matrix<double, 5, 3> _fcA;
            _fcA.setZero();
            _fcA << 1, 0, u, -1, 0, u, 0, 1, u, 0, -1, u, 0, 0, 1;
            cA.block<5, 3>(0, 0) = _fcA;
            cA.block<5, 3>(5, 3) = _fcA;
            cA.block<5, 3>(10, 6) = _fcA;
            cA.block<5, 3>(15, 9) = _fcA;
            Eigen::Matrix<double,8,16> _tcA;
            _tcA.setZero();
            for (int i = 0; i < 4; i++)
            {
                _tcA(0 + 2 * i, 2 + 3 * i) = u;
                _tcA(1 + 2 * i, 2 + 3 * i) = u;
                _tcA(0 + 2 * i, 12 + i) = 1. / bodyObject->legs[i]->Reff;
                _tcA(1 + 2 * i, 12 + i) = - 1. / bodyObject->legs[i]->Reff;
            }
            Alb.setZero();
            Aub.setConstant(100000.);
        }

        // 更新机器人动力学方程（描述为 left*[f] = right）
        void updateDynamic() override;
        void updateDynamic(const Eigen::Matrix3d& _slope);
        // 导入权重参数
        void importWeight(const VectorXd& _Q, const VectorXd& _F, const VectorXd& _R, const VectorXd& _W) override;
        // 导入力参数
        void importForce(const double& _force_cz, const double& _force_cxy, const double& _u, const double& _tau_c);
        // 更新当前状态
        void updateBalanceState() override;
        // 设置终端目标
        void setPositionTarget(const Vector3d& _p, const Vector3d& _r, const Vector4d& _pe);
        void setVelocityTarget(const Vector3d& _p_dot, const Vector3d& _r_dot, const Vector4d& _pe_dot);
        // 执行mpc控制器
        void mpc_adjust(const VectorX<bool>& _enList) override;
        // 设置接触约束
        void setContactConstrain(const Vector4i& _contact,const Eigen::Matrix<double, 3, 4>& _swingForce);
        void setContactConstrain(const Vector4i& _contact, const Eigen::Matrix<double, 3, 4>& _swingForce, const Eigen::Matrix3d& _slope);
        // 直接设置四足输出
        void setLegsForce(const Eigen::Matrix<double, 3, 4>& _force, const Eigen::Vector4d& _tau);
        // 四足腾空处理
        void contactDeal(const VectorXd& _oriQ, const double _ffRatio, const double _stRatio);
    };

    void QpwPVCtrl::updateDynamic()
    {
        for (int i = 0; i < 4; i++)
        {
            // 四足接触点位置反对称矩阵的计算
            //dynamicLeft.block<3, 3>(0, i * 3) = Rsb_c;
            Vector3d Pbi = Vector3d::Zero();
            Pbi = bodyObject->Rsb_c * (bodyObject->currentBodyState.leg_b[i].Position - bodyObject->P);
            //Pbi = (currentBodyState.leg_b[i].Position - Pb);
            dynamicLeft.block<3, 3>(0, i * 3) = fftauRatio[i] * Eigen::Matrix3d::Identity();
            dynamicLeft.block<3, 3>(3, i * 3) = fftauRatio[i] * bodyObject->v3_to_m3(Pbi);

            // 轮相关物理参数更新
            dynamicRight(6 + i, 6 + i) = bodyObject->legs[i]->Ic[3](1) / pow(bodyObject->legs[i]->Reff, 2) + bodyObject->legs[i]->Mc[3];
            dynamicLeft(6 + i, 12 + i) = 1. / bodyObject->legs[i]->Reff;
        }
        dynamicRight.block<3, 3>(0, 0) = bodyObject->M;
        dynamicRight.block<3, 3>(3, 3) = bodyObject->Rsb_c * bodyObject->I * bodyObject->Rsb_c.transpose()/* * bodyObject->dEuler2W*/;

        // 为了减少计算量，轮速规划是基于车体坐标系的，因此要对世界坐标系的力做映射
        Eigen::Matrix<double, 4, 3> s = Eigen::Matrix<double, 4, 3>::Zero();
        //s(0, 0) = -fftauRatio[0];
        s(0, 0) = -1;
        dynamicLeft.block(6, 0, 4, 3) = s * bodyObject->Rsbh_c.transpose();
        s(0, 0) = 0;
        /*s(1, 0) = -fftauRatio[1];*/
        s(1, 0) = -1;
        dynamicLeft.block(6, 3, 4, 3) = s * bodyObject->Rsbh_c.transpose();
        s(1, 0) = 0;
        //s(2, 0) = -fftauRatio[2];
        s(2, 0) = -1;
        dynamicLeft.block(6, 6, 4, 3) = s * bodyObject->Rsbh_c.transpose();
        s(2, 0) = 0;
        /*s(3, 0) = -fftauRatio[3];*/
        s(3, 0) = -1;
        dynamicLeft.block(6, 9, 4, 3) = s * bodyObject->Rsbh_c.transpose();
    }

    void QpwPVCtrl::updateDynamic(const Eigen::Matrix3d& _slope)
    {
        for (int i = 0; i < 4; i++)
        {
            // 四足接触点位置反对称矩阵的计算
            //dynamicLeft.block<3, 3>(0, i * 3) = Rsb_c;
            Vector3d Pbi = Vector3d::Zero();
            Pbi = bodyObject->Rsb_c * (bodyObject->currentBodyState.leg_b[i].Position - bodyObject->P);
            //Pbi = (currentBodyState.leg_b[i].Position - Pb);
            dynamicLeft.block<3, 3>(0, i * 3) = fftauRatio[i] * Eigen::Matrix3d::Identity();
            dynamicLeft.block<3, 3>(3, i * 3) = fftauRatio[i] * bodyObject->v3_to_m3(Pbi);

            // 轮相关物理参数更新
            dynamicRight(6 + i, 6 + i) = bodyObject->legs[i]->Ic[3](1) / pow(bodyObject->legs[i]->Reff, 2) + bodyObject->legs[i]->Mc[3];
            dynamicLeft(6 + i, 12 + i) = 1. / bodyObject->legs[i]->Reff;
        }
        dynamicRight.block<3, 3>(0, 0) = bodyObject->M;
        dynamicRight.block<3, 3>(3, 3) = bodyObject->Rsb_c * bodyObject->I * bodyObject->Rsb_c.transpose()/* * bodyObject->dEuler2W*/;

        // 为了减少计算量，轮速规划是基于车体坐标系的，因此要对世界坐标系的力做映射
        Eigen::Matrix<double, 4, 3> s = Eigen::Matrix<double, 4, 3>::Zero();
        //s(0, 0) = -fftauRatio[0];
        s(0, 0) = -1;
        dynamicLeft.block(6, 0, 4, 3) = s * bodyObject->Rsbh_c.transpose() * _slope.transpose();
        s(0, 0) = 0;
        /*s(1, 0) = -fftauRatio[1];*/
        s(1, 0) = -1;
        dynamicLeft.block(6, 3, 4, 3) = s * bodyObject->Rsbh_c.transpose() * _slope.transpose();
        s(1, 0) = 0;
        //s(2, 0) = -fftauRatio[2];
        s(2, 0) = -1;
        dynamicLeft.block(6, 6, 4, 3) = s * bodyObject->Rsbh_c.transpose() * _slope.transpose();
        s(2, 0) = 0;
        /*s(3, 0) = -fftauRatio[3];*/
        s(3, 0) = -1;
        dynamicLeft.block(6, 9, 4, 3) = s * bodyObject->Rsbh_c.transpose() * _slope.transpose();
    }

    void QpwPVCtrl::importWeight(const VectorXd& _Q, const VectorXd& _F, const VectorXd& _R, const VectorXd& _W)
    {
        for (int i = 0; i < 23; i++)
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

    void QpwPVCtrl::importForce(const double& _force_cz, const double& _force_cxy, const double& _u, const double& _tau_c)
    {
        this->force_cz = _force_cz;
        this->force_cxy = _force_cxy;
        this->u = _u;
        this->tau_c = _tau_c;
    }

    void QpwPVCtrl::updateBalanceState()
    {
        static multiCircle angC = multiCircle(3.1415926);
        currentBalanceState.p = bodyObject->currentWorldState.dist;
        currentBalanceState.p_dot = bodyObject->currentWorldState.linVel_xyz;
        currentBalanceState.r(0) = bodyObject->currentBodyState.Ang_xyz(0);
        currentBalanceState.r(1) = bodyObject->currentBodyState.Ang_xyz(1);
        currentBalanceState.r(2) = angC.f(bodyObject->currentBodyState.Ang_xyz(2));
        /*currentBalanceState.r_dot = bodyObject->Rsb_c * bodyObject->dEuler2W.inverse() * bodyObject->currentBodyState.angVel_xyz;*/
        currentBalanceState.r_dot = bodyObject->currentWorldState.angVel_xyz;
        for (int i = 0; i < 4; i++)
        {
            currentBalanceState.pe(i) = (bodyObject->Rsbh_c.transpose() * bodyObject->Rsb_c * bodyObject->currentBodyState.leg_b[i].Position)(0);
            currentBalanceState.pe_dot(i) = (bodyObject->Rsbh_c.transpose() * bodyObject->currentWorldState.leg_s[i].VelocityW)(0);
        }
        //std::cout << "cbpe: \n" << currentBalanceState.pe_dot << std::endl;
    }

    void QpwPVCtrl::setPositionTarget(const Vector3d& _p, const Vector3d& _r, const Vector4d& _pe)
    {
        targetBalanceState.p = _p;
        targetBalanceState.r = _r;
        targetBalanceState.pe = _pe;
    }

    void QpwPVCtrl::setVelocityTarget(const Vector3d& _p_dot, const Vector3d& _r_dot, const Vector4d& _pe_dot)
    {
        targetBalanceState.p_dot = _p_dot;
        targetBalanceState.r_dot = _r_dot;
        targetBalanceState.pe_dot = _pe_dot;
    }

    void QpwPVCtrl::mpc_adjust(const VectorX<bool>& _enList)
    {
        A.block(0, 10, 3, 3) = Eigen::Matrix3d::Identity();
        A.block(3, 13, 3, 3) = bodyObject->Rsb_c.transpose();
        A.block(6, 16, 4, 4) = Eigen::Matrix4d::Identity();
        A.block(10, 20, 3, 3) = Eigen::Matrix3d::Identity();
        B.block(10,0,10,16) = dynamicRight.inverse() * dynamicLeft;
        balanceController.setConstrain(lb, ub);
        balanceController.setBoxConstrain(cA, Alb, Aub);
        y.block<3, 1>(0, 0) = targetBalanceState.p;
        y.block<3, 1>(3, 0) = targetBalanceState.r;
        y.block<4, 1>(6, 0) = targetBalanceState.pe;
        y.block<3, 1>(10, 0) = targetBalanceState.p_dot;
        y.block<3, 1>(13, 0) = targetBalanceState.r_dot;
        y.block<4, 1>(16, 0) = targetBalanceState.pe_dot;
        y.block<3, 1>(20, 0) = g;
        x.block<3, 1>(0, 0) = currentBalanceState.p;
        x.block<3, 1>(3, 0) = currentBalanceState.r;
        x.block<4, 1>(6, 0) = currentBalanceState.pe;
        x.block<3, 1>(10, 0) = currentBalanceState.p_dot;
        x.block<3, 1>(13, 0) = currentBalanceState.r_dot;
        x.block<4, 1>(16, 0) = currentBalanceState.pe_dot;
        x.block<3, 1>(20, 0) = g;
        balanceController.mpc_update(y, x, 100, 0.02);
        balanceController.mpc_init(A, B, Q, F, R, W, dt);
        balanceController.mpc_solve(0);
        for (int i = 0; i < 4; i++)
        {
            this->mpcOut.block(0, i, 3, 1) = -bodyObject->Rsb_c.transpose() * balanceController.getOutput().block<3, 1>(3 * i, 0);
            this->mpcOut(3, i) = balanceController.getOutput()(12 + i, 0);
        }
    }

    void QpwPVCtrl::setContactConstrain(const Vector4i& _contact,const Eigen::Matrix<double, 3, 4>& _swingForce)
    {
        Eigen::Matrix<double, 5, 3> _fcA;
        _fcA.setZero();
        Eigen::Matrix<double, 5, 1> _Aub;
        _Aub.setZero();
        Eigen::Matrix<double, 8, 16> _tcA;
        _tcA.setZero();
        // 若该腿不触地，则清零约束矩阵
        for (int i = 0; i < 4; i++)
        {
            if (_contact(i) == 1)
            {
                _fcA << 1, 0, u, -1, 0, u, 0, 1, u, 0, -1, u, 0, 0, 1;
                _Aub.setConstant(100000.);
                _tcA(0 + 2 * i, 2 + 3 * i) = u;
                _tcA(1 + 2 * i, 2 + 3 * i) = u;
                _tcA(0 + 2 * i, 12 + i) = 1. / bodyObject->legs[i]->Reff;
                _tcA(1 + 2 * i, 12 + i) = -1. / bodyObject->legs[i]->Reff;
            }
            else
            {
                _fcA << 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0;
                _Aub.setZero();
                _tcA(0 + 2 * i, 2 + 3 * i) = 0;
                _tcA(1 + 2 * i, 2 + 3 * i) = 0;
                _tcA(0 + 2 * i, 12 + i) = 1.;
                _tcA(1 + 2 * i, 12 + i) = 0;
                // 如果不接触那么就更新摆动力到mpc控制器中
                this->balanceController.updateLastU(_swingForce.col(i), i * 3, 3);
            }
            this->cA.block<5, 3>(5 * i, 3 * i) = _fcA;
            this->Aub.block<5, 1>(5 * i, 0) = _Aub;
        }
        this->cA.block(20, 0, 8, 16) = _tcA;
        for (int i = 0; i < 4; i++)
        {
            lb(3 * i, 0) = -force_cxy;
            lb(3 * i + 1, 0) = -force_cxy;
            lb(3 * i + 2, 0) = -force_cz;
            ub(3 * i, 0) = force_cxy;
            ub(3 * i + 1, 0) = force_cxy;
            ub(3 * i + 2, 0) = force_cz;
        }
        lb.block(12, 0, 4, 1).setConstant(-tau_c);
        ub.block(12, 0, 4, 1).setConstant(tau_c);
    }

    void QpwPVCtrl::setContactConstrain(const Vector4i& _contact, const Eigen::Matrix<double, 3, 4>& _swingForce, const Eigen::Matrix3d& _slope)
    {
        Eigen::Matrix<double, 5, 3> _fcA;
        _fcA.setZero();
        Eigen::Matrix<double, 5, 1> _Aub;
        _Aub.setZero();
        Eigen::Matrix<double, 8, 16> _tcA;
        _tcA.setZero();
        // 若该腿不触地，则清零约束矩阵
        for (int i = 0; i < 4; i++)
        {
            if (_contact(i) == 1)
            {
                _fcA << 1, 0, u, -1, 0, u, 0, 1, u, 0, -1, u, 0, 0, 1;
                _fcA = _fcA * _slope.transpose();
                _Aub.setConstant(100000.);
                _tcA(0 + 2 * i, 2 + 3 * i) = u;
                _tcA(1 + 2 * i, 2 + 3 * i) = u;
                _tcA.block(2 * i, 3 * i, 2, 3) = _tcA.block(2 * i, 3 * i, 2, 3) * bodyObject->Rsbh_c.transpose() * _slope.transpose();
                _tcA(0 + 2 * i, 12 + i) = 1. / bodyObject->legs[i]->Reff;
                _tcA(1 + 2 * i, 12 + i) = -1. / bodyObject->legs[i]->Reff;
                /*_tcA(0 + 2 * i, 2 + 3 * i) = 0;
                _tcA(1 + 2 * i, 2 + 3 * i) = 0;
                _tcA(0 + 2 * i, 12 + i) = 0.;
                _tcA(1 + 2 * i, 12 + i) = 0;*/
            }
            else
            {
                _fcA << 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0;
                _Aub.setZero();
                _tcA(0 + 2 * i, 2 + 3 * i) = 0;
                _tcA(1 + 2 * i, 2 + 3 * i) = 0;
                _tcA(0 + 2 * i, 12 + i) = 1.;
                _tcA(1 + 2 * i, 12 + i) = 0;
                // 如果不接触那么就更新摆动力到mpc控制器中
                this->balanceController.updateLastU(_swingForce.col(i), i * 3, 3);
            }
            this->cA.block<5, 3>(5 * i, 3 * i) = _fcA;
            this->Aub.block<5, 1>(5 * i, 0) = _Aub;
        }
        this->cA.block(20, 0, 8, 16) = _tcA;
        for (int i = 0; i < 4; i++)
        {
            lb(3 * i, 0) = -force_cxy;
            lb(3 * i + 1, 0) = -force_cxy;
            lb(3 * i + 2, 0) = -force_cz;
            ub(3 * i, 0) = force_cxy;
            ub(3 * i + 1, 0) = force_cxy;
            ub(3 * i + 2, 0) = force_cz;
        }
        lb.block(12, 0, 4, 1).setConstant(-tau_c);
        ub.block(12, 0, 4, 1).setConstant(tau_c);
    }

    void QpwPVCtrl::setLegsForce(const Eigen::Matrix<double, 3, 4>& _force, const Eigen::Vector4d& _tau)
    {
        for (int i = 0; i < 4; i++)
        {
            this->legsCtrl[i]->setEndForceTar(_force.col(i));
            this->legsCtrl[i]->setEndTauTar(_tau(i));
        }
    }

    void QpwPVCtrl::contactDeal(const VectorXd& _oriQ, const double _ffRatio, const double _stRatio)
    {
        this->Q = _oriQ.asDiagonal();
        if (_stRatio < 1.)
        {
            for (int i = 0; i < 4; i++)
            {
                this->Q(6 + i, 6 + i) = _oriQ(6 + i) * (1e-10 + bodyObject->est->_ctTrust[i]);
                this->Q(16 + i, 16 + i) = _oriQ(16 + i) * (1 + 1e1 * (1 - bodyObject->est->_ctTrust[i]));
                this->fftauRatio[i] = _ffRatio * bodyObject->est->_ctTrust[i];
                /*this->Q(6 + i, 6 + i) = _oriQ(6 + i) * (1e-10 + bodyObject->mixContact(i));
                this->Q(16 + i, 16 + i) = _oriQ(16 + i) * (1 + 1e1 * (1 - bodyObject->mixContact(i)));*/
                //this->fftauRatio[i] = _ffRatio * bodyObject->mixContact(i);
            }
        }
        else
        {
            for (int i = 0; i < 4; i++)
            {
                this->fftauRatio[i] = 1;
            }
        }
        this->F = this->Q;
    }
}
