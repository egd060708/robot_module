#pragma once

#include <Eigen/Dense>
#include "mpcCal.h"
#include <iostream>

using namespace Eigen;

namespace Quadruped {

	class PlanBase {
	protected:
		/* 机身的轨迹规划 */
		Vector3d Pbt,Pbc,Vbt,Vbc;
		/* 足端的轨迹规划 */
		Eigen::Matrix<double, 3, 4> Plt, Plc, Vlt, Vlc;
		/* 规划结果 */
		Vector3d Pbr, Vbr;
		Eigen::Matrix<double, 3, 4> Plr, Vlr;
	public:
		void updateBodyHighLevelTar(const Vector3d& _Pbt, const Vector3d& _Vbt)
		{
			this->Pbt = _Pbt;
			this->Vbt = _Vbt;
		}
		void updateFootHighLevelTar(const Eigen::Matrix<double, 3, 4>& _Plt, const Eigen::Matrix<double, 3, 4>& _Vlt)
		{
			this->Plt = _Plt;
			this->Vlt = _Vlt;
		}
		void updateBodyState(const Vector3d& _Pbc, const Vector3d& _Vbc)
		{
			this->Pbc = _Pbc;
			this->Vbc = _Vbc;
		}
		void updateFootState(const Eigen::Matrix<double, 3, 4>& _Plc, const Eigen::Matrix<double, 3, 4>& _Vlc)
		{
			this->Plc = _Plc;
			this->Vlc = _Vlc;
		}
		virtual void useBodyPlan() = 0;
		virtual void useFootPlan() = 0;
		const Vector3d& getBodyPlanPosition()
		{
			return this->Pbr;
		}
		const Vector3d& getBodyPlanVelocity()
		{
			return this->Vbr;
		}
		const Eigen::Matrix<double, 3, 4>& getFootPlanPosition()
		{
			return this->Plr;
		}
		const Eigen::Matrix<double, 3, 4>& getFootPlanVelocity()
		{
			return this->Vlr;
		}
		void setInitBodyPlanPosition(const Eigen::Vector3d& _Pbr)
		{
			this->Pbr = _Pbr;
		}
		void setInitBodyPlanVelocity(const Eigen::Vector3d& _Vbr)
		{
			this->Vbr = _Vbr;
		}
	};

	class BalanceBodyPlanning :public PlanBase
	{
	private:
		/* 翻滚轴参数 */
		Vector4d tumblingAxis;
		/* 翻滚点参数 */
		Eigen::Matrix<double, 3, 4> tumblingPoint;
		/* 平衡锥约束系数，通过自定义稳定裕度以及支撑多边形找到保守的平衡追约束 */
		double balanceMargin, balanceCone;
		/* 轨迹规划mpc生成器 */
		mpcCal<6, 3, 4, 1, 10> traGenerator;
		/* 预测时间步长 */
		double step_t = 0;
		/* mpc输出观测 */
		Eigen::Vector3d acc;

		/* 获取翻滚轴以及翻滚点 */
		void getTumbling()
		{
			for (int i = 0; i < 4; i++)
			{
				if (i == 3)
				{
					double temp = -(this->Plc.col(0) - this->Plc.col(i)).transpose() * (this->Plc.col(i) - this->Pbc);
					tumblingAxis(i) = temp / (this->Plc.col(0) - this->Plc.col(i)).squaredNorm();
					tumblingPoint.col(i) = this->Plc.col(i) + tumblingAxis(i) * (this->Plc.col(0) - this->Plc.col(i));
				}
				else
				{
					double temp = -(this->Plc.col(i + 1) - this->Plc.col(i)).transpose() * (this->Plc.col(i) - this->Pbc);
					tumblingAxis(i) = temp / (this->Plc.col(i + 1) - this->Plc.col(i)).squaredNorm();
					tumblingPoint.col(i) = this->Plc.col(i) + tumblingAxis(i) * (this->Plc.col(i + 1) - this->Plc.col(i));
				}
			}
			/*std::cout << "Axis: " << std::endl;
			std::cout << tumblingAxis << std::endl;
			std::cout << "Point: " << std::endl;
			std::cout << tumblingPoint << std::endl;*/
		}

		/* 找到平衡锥约束参数 */
		void findBalanceCone()
		{
			Eigen::Vector3d g(0, 0, -9.81);
			// 找到最严格的约束点(计算最小夹角)
			double temp = g.dot(this->tumblingPoint.col(0) - this->Pbc) / (g.norm() + (this->tumblingPoint.col(0) - this->Pbc).norm());
			int min_idx = 0;
			for (int i = 0; i < 4; i++)
			{
				double temp1 = g.dot(this->tumblingPoint.col(i) - this->Pbc) / (g.norm() + (this->tumblingPoint.col(i) - this->Pbc).norm());
				if (temp1 < temp)
				{
					min_idx = i;
					temp = temp1;
				}
			}
			balanceCone = balanceMargin * (this->tumblingPoint.col(min_idx) - this->Pbc).block(0, 0, 2, 1).norm() / abs((this->tumblingPoint.col(min_idx) - this->Pbc)(2));
			//std::cout << balanceCone << std::endl;
		}

		/* 更新不等式约束 */
		void updateIeqConstrain()
		{
			// 平衡约束转化为线性锥约束
			Eigen::Matrix<double, 4, 3> cA;
			cA.setZero();
			cA << 1, 0, this->balanceCone, -1, 0, this->balanceCone, 0, 1, this->balanceCone, 0, -1, this->balanceCone;
			Eigen::Vector4d Aub, Alb;
			Aub.setConstant(-(this->balanceCone * -9.81));
			Alb.setConstant(-100000.);
			this->traGenerator.setBoxConstrain(cA, Alb, Aub);
			/*std::cout << "cA: " << std::endl;
			std::cout << cA << std::endl;
			std::cout << "Alb: " << std::endl;
			std::cout << Alb << std::endl;
			std::cout << "Aub: " << std::endl;
			std::cout << Aub << std::endl;*/
		}

		/* 求解优化问题并得到输出轨迹 */
		void getSolution()
		{
			Eigen::Vector<double, 6> y, x;
			Eigen::Vector3d Py(this->Pbr(0), this->Pbr(1), this->Pbt(2));
			x.block(0, 0, 3, 1) = this->Pbr;
			x.block(3, 0, 3, 1) = this->Vbr;
			y.block(0, 0, 3, 1) = Py;
			y.block(3, 0, 3, 1) = this->Vbt;
			/*std::cout << "x: " << std::endl;
			std::cout << x << std::endl;
			std::cout << "y: " << std::endl;
			std::cout << y << std::endl;*/
			this->traGenerator.mpc_update(y, x, 100, 1.);
			this->traGenerator.mpc_solve();
			this->acc = this->traGenerator.getOutput();
			this->Vbr = this->Vbr + 0.002 * this->traGenerator.getOutput();
			this->Pbr = this->Pbr + 0.002 * this->Vbr;
		}

	public:
		/* 构造函数 */
		BalanceBodyPlanning():traGenerator(PL_LOW) {
			// 初始化参数
			this->balanceCone = 0.5;
			this->balanceMargin = 0.5;
			this->Pbr.setZero();
			this->Vbr.setZero();
		}
		/* 使用机身规划器 */
		void useBodyPlan() override
		{
			this->getTumbling();
			this->findBalanceCone();
			this->updateIeqConstrain();
			this->getSolution();
		}
		/* 使用足端规划器 */
		void useFootPlan() override
		{

		}
		/* 更新规划器参数 */
		void updateParams(const double& _t, const double& _balanceMargin, const Eigen::Matrix<double, 6, 6>& _Q, const Eigen::Matrix<double, 6, 6>& _F, const Eigen::Matrix<double, 3, 3>& _R, const Eigen::Matrix<double,3,3>& _W)
		{
			this->step_t = _t;
			this->balanceMargin = _balanceMargin;
			// 初始化轨迹生成器
			MatrixXd A, B;
			A.resize(6, 6);
			B.resize(6, 3);
			A.setZero();
			B.setZero();
			A.block(0, 0, 3, 3).setIdentity();
			Eigen::Vector3d st;
			st.setConstant(_t);
			A.block(0, 3, 3, 3) = st.asDiagonal();
			A.block(3, 3, 3, 3).setIdentity();
			B.block(3, 0, 3, 3) = st.asDiagonal();
			/*std::cout << "A: " << std::endl;
			std::cout << A << std::endl;
			std::cout << "B: " << std::endl;
			std::cout << B << std::endl;*/
			traGenerator.mpc_init(A, B, _Q, _F, _R, _W, 0);
		}
		/* 更新等式约束 */
		void setEqConstrain(const Vector3d& _aMax)
		{
			// 加速度上下限(对向量每个元素取绝对值操作)
			Eigen::Vector3d ub = _aMax.cwiseAbs();
			Eigen::Vector3d lb = -_aMax.cwiseAbs();
			traGenerator.setConstrain(lb, ub);
		}
		/* 获取mpc原始输出 */
		const Eigen::Vector3d& getMpcOriOut()
		{
			return this->acc;
		}
	};

}