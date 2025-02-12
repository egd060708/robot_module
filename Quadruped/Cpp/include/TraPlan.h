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
		PlanBase()
		{
			this->Pbt.setZero();
			this->Pbc.setZero();
			this->Vbt.setZero();
			this->Vbc.setZero();
			this->Plt.setZero();
			this->Plc.setZero();
			this->Vlt.setZero();
			this->Vlc.setZero();
			this->Pbr.setZero();
			this->Vbr.setZero();
			this->Plr.setZero();
			this->Vlr.setZero();
		}
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
		void setInitFootPlanPosition(const Eigen::Matrix<double, 3, 4>& _Plr)
		{
			this->Plr = _Plr;
		}
		void setInitFootPlanVelocity(const Eigen::Matrix<double, 3, 4>& _Vlr)
		{
			this->Vlr = _Vlr;
		}
	};

	class SupportPolygon :public PlanBase
	{
	protected:
		/* 翻滚轴参数 */
		Vector4d tumblingAxis;
		/* 翻滚点参数 */
		Eigen::Matrix<double, 3, 4> tumblingPoint;
		/* 平衡锥约束系数，通过自定义稳定裕度以及支撑多边形找到保守的平衡追约束 */
		double balanceMargin, balanceCone;
		double balanceMarginX, balanceMarginY;
		double balanceConeX, balanceConeY;

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
		}

		/* 找到平衡锥约束参数 */
		void findBalanceCone()
		{
			Eigen::Vector3d g(0, 0, -9.81);
			// 找到最严格的约束点(计算最小夹角)，以及相对于x轴和y轴的最严格约束点
			double temp = g.dot(this->tumblingPoint.col(0) - this->Pbc) / (g.norm() + (this->tumblingPoint.col(0) - this->Pbc).norm());
			double tempx = temp;
			double tempy = g.dot(this->tumblingPoint.col(1) - this->Pbc) / (g.norm() + (this->tumblingPoint.col(1) - this->Pbc).norm());
			int min_idx = 0, min_idxX = 0, min_idxY = 1;
			for (int i = 0; i < 4; i++)
			{
				double temp1 = g.dot(this->tumblingPoint.col(i) - this->Pbc) / (g.norm() + (this->tumblingPoint.col(i) - this->Pbc).norm());
				if (temp1 < temp)
				{
					min_idx = i;
					temp = temp1;
				}
				if (i == 0 || i == 2)
				{
					if (temp1 < tempx)
					{
						min_idxX = i;
						tempx = temp1;
					}
				}
				else if (i == 1 || i == 3)
				{
					if (temp1 < tempy)
					{
						min_idxY = i;
						tempy = temp1;
					}
				}
			}
			balanceCone = balanceMargin * (this->tumblingPoint.col(min_idx) - this->Pbc).block(0, 0, 2, 1).norm() / abs((this->tumblingPoint.col(min_idx) - this->Pbc)(2));
			balanceConeX = balanceMarginX * (this->tumblingPoint.col(min_idxX) - this->Pbc).block(0, 0, 2, 1).norm() / abs((this->tumblingPoint.col(min_idxX) - this->Pbc)(2));
			balanceConeY = balanceMarginY * (this->tumblingPoint.col(min_idxY) - this->Pbc).block(0, 0, 2, 1).norm() / abs((this->tumblingPoint.col(min_idxY) - this->Pbc)(2));
		}

	public:
		/* 构造函数 */
		SupportPolygon(){
			// 初始化参数
			this->balanceCone = 0.5;
			this->balanceMargin = 0.5;
			this->balanceMarginX = 0.5;
			this->balanceMarginY = 0.5;

		}
		/* 更新支撑多边形参数 */
		void updateSupportPolygon()
		{
			this->getTumbling();
			this->findBalanceCone();
		}
	};

	class BalanceDividePlanning : public SupportPolygon
	{
	private:
		/* 针对机身的轨迹规划mpc生成器 */
		mpcCal<6, 3, 16, 1, 10> traBodyGenerator;
		/* 预测时间步长 */
		double body_step_t = 0;
		double foot_step_t = 0;
		/* 机身mpc输出观测 */
		Eigen::Vector3d bodyAcc;
		/* 针对足端的规划生成器 */
		mpcCal<16, 8, 16, 1, 10> traFootGenerator;
		/* 足端mpc输出观测 */
		Eigen::Matrix<double, 3, 4> footAcc;
		/* 坐标转换 */
		Eigen::Matrix3d Rsb_c;

		/* 更新机身不等式约束 */
		void updateBodyIeqConstrain(const Vector3d& _aMax)
		{
			// 平衡约束转化为线性锥约束
			Eigen::Matrix<double, 16, 3> cA;
			cA.setZero();
			cA << 1, 0, this->balanceCone, -1, 0, this->balanceCone, 0, 1, this->balanceCone, 0, -1, this->balanceCone, \
				1, 1, 1.414 * this->balanceCone, -1, 1, 1.414 * this->balanceCone, -1, -1, 1.414 * this->balanceCone, 1, -1, 1.414 * this->balanceCone, \
				1, 1, -1.414 * this->balanceCone, -1, 1, -1.414 * this->balanceCone, -1, -1, -1.414 * this->balanceCone, 1, -1, 1.414 * this->balanceCone, \
				1, 1, 0, 1, -1, 0, -1, 1, 0, -1, -1, 0;
			Eigen::Vector<double, 16> Aub, Alb;
			Aub.block(0, 0, 4, 1).setConstant(-(this->balanceCone * -9.81));
			Alb.block(0, 0, 4, 1).setConstant(-100000.);
			Aub.block(4, 0, 4, 1).setConstant(-1.414 * this->balanceCone * -9.81);
			Alb.block(4, 0, 4, 1).setConstant(-100000.);
			Aub.block(8, 0, 4, 1).setConstant(100000.);
			Alb.block(8, 0, 4, 1).setConstant(1.414 * this->balanceCone * -9.81);
			Aub.block(12, 0, 4, 1).setConstant(_aMax.segment(0, 2).norm());
			Alb.block(12, 0, 4, 1).setConstant(-_aMax.segment(0, 2).norm());
			this->traBodyGenerator.setBoxConstrain(cA, Alb, Aub);
		}

		/* 求解机身轨迹优化问题并得到输出轨迹 */
		void getBodySolution()
		{
			static Eigen::Vector3d last_Vbr = this->Vbr;
			Eigen::Vector<double, 6> y, x;
			Eigen::Vector3d Py(this->Pbr(0), this->Pbr(1), this->Pbt(2));
			x.block(0, 0, 3, 1) = this->Pbr;
			x.block(3, 0, 3, 1) = this->Vbr;
			y.block(0, 0, 3, 1) = this->Pbt;
			y.block(3, 0, 3, 1) = this->Vbt;
			this->traBodyGenerator.mpc_update(y, x, 100, 1.);
			this->traBodyGenerator.mpc_solve();
			this->bodyAcc = this->traBodyGenerator.getOutput();
			this->Vbr = this->Vbr + 0.002 * this->traBodyGenerator.getOutput();
			this->Pbr = this->Pbr + 0.5 * 0.002 * (last_Vbr + this->Vbr);
			last_Vbr = this->Vbr;
		}

		/* 更新足端不等式约束 */
		void updateFootIeqConstrain()
		{
			// 对足端的加速度边界不等式约束设置所有边界
			double edgeX1, edgeX2, edgeY1, edgeY2;
			Eigen::Vector3d realBodyAcc = this->Rsb_c.transpose() * this->bodyAcc;
			edgeX1 = realBodyAcc(0) + this->balanceConeX * (realBodyAcc(2) - 9.81);
			edgeX2 = -realBodyAcc(0) + this->balanceConeX * (realBodyAcc(2) - 9.81);
			edgeY1 = realBodyAcc(1) + this->balanceConeY * (realBodyAcc(2) - 9.81);
			edgeY2 = -realBodyAcc(1) + this->balanceConeY * (realBodyAcc(2) - 9.81);
			Eigen::Matrix<double, 16, 8> cA;
			cA.setZero();
			cA(0, 0) = 1; cA(1, 0) = 1;
			cA(2, 1) = 1; cA(3, 1) = 1;
			cA(4, 2) = 1; cA(5, 2) = 1;
			cA(6, 3) = -1; cA(7, 3) = -1;
			cA(8, 4) = -1; cA(9, 4) = -1;
			cA(10, 5) = 1; cA(11, 5) = 1;
			cA(12, 6) = -1; cA(13, 6) = -1;
			cA(14, 7) = -1; cA(15, 7) = -1;
			Eigen::Vector<double, 16> Aub, Alb;
			Aub.setConstant(100000.);
			for (int i = 0; i < 4; i++)
			{
				Alb(4 * i) = edgeX1;
				Alb(4 * i + 1) = edgeX2;
				Alb(4 * i + 2) = edgeY1;
				Alb(4 * i + 3) = edgeY2;
			}
			this->traFootGenerator.setBoxConstrain(cA, Alb, Aub);
		}

		/* 求解足端轨迹优化问题并得到输出轨迹 */
		void getFootSolution()
		{
			static Eigen::Matrix<double, 3, 4> last_Vlr = this->Vlr;
			Eigen::Vector<double, 16> y, x;
			x.block(0, 0, 2, 1) = this->Plr.col(0).segment(0, 2);
			x.block(2, 0, 2, 1) = this->Plr.col(1).segment(0, 2);
			x.block(4, 0, 2, 1) = this->Plr.col(2).segment(0, 2);
			x.block(6, 0, 2, 1) = this->Plr.col(3).segment(0, 2);
			x.block(8, 0, 2, 1) = this->Vlr.col(0).segment(0, 2);
			x.block(10, 0, 2, 1) = this->Vlr.col(1).segment(0, 2);
			x.block(12, 0, 2, 1) = this->Vlr.col(2).segment(0, 2);
			x.block(14, 0, 2, 1) = this->Vlr.col(3).segment(0, 2);
			y.block(0, 0, 2, 1) = this->Plt.col(0).segment(0, 2);
			y.block(2, 0, 2, 1) = this->Plt.col(1).segment(0, 2);
			y.block(4, 0, 2, 1) = this->Plt.col(2).segment(0, 2);
			y.block(6, 0, 2, 1) = this->Plt.col(3).segment(0, 2);
			y.block(8, 0, 8, 1).setZero();
			this->traFootGenerator.mpc_update(y, x, 100, 1.);
			this->traFootGenerator.mpc_solve();
			this->footAcc.setZero();
			this->footAcc.block(0, 0, 2, 1) = this->traFootGenerator.getOutput().block(0, 0, 2, 1);
			this->footAcc.block(0, 1, 2, 1) = this->traFootGenerator.getOutput().block(2, 0, 2, 1);
			this->footAcc.block(0, 2, 2, 1) = this->traFootGenerator.getOutput().block(4, 0, 2, 1);
			this->footAcc.block(0, 3, 2, 1) = this->traFootGenerator.getOutput().block(6, 0, 2, 1);
			this->Vlr = this->Vlr + 0.002 * this->footAcc;
			this->Plr = this->Plr + 0.5 * 0.002 * (last_Vlr + this->Vlr);
			last_Vlr = this->Vlr;
		}
	public:
		/* 构造函数 */
		BalanceDividePlanning():traBodyGenerator(PL_NONE),traFootGenerator(PL_NONE)  {
			// 初始化参数
			this->Pbr.setZero();
			this->Vbr.setZero();
			this->Plr.setZero();
			this->Vlr.setZero();
			this->Rsb_c.setIdentity();
		}
		/* 规划器预热 */
		void warmUp(const Eigen::Matrix3d& _Rsbc)
		{
			this->Rsb_c = _Rsbc;
			this->updateSupportPolygon();
		}
		/* 使用机身规划器 */
		void useBodyPlan()
		{
			//this->updateBodyIeqConstrain();
			this->getBodySolution();
		}
		/* 使用足端规划器 */
		void useFootPlan()
		{
			this->updateFootIeqConstrain();
			this->getFootSolution();
		}
		/* 更新规划器参数 */
		void updateBodyParams(const double& _t, const double& _balanceMargin, const Eigen::Matrix<double, 6, 6>& _Q, const Eigen::Matrix<double, 6, 6>& _F, const Eigen::Matrix<double, 3, 3>& _R, const Eigen::Matrix<double, 3, 3>& _W)
		{
			this->body_step_t = _t;
			this->balanceMargin = _balanceMargin;
			// 初始化轨迹生成器
			MatrixXd A, B;
			A.resize(6, 6);
			B.resize(6, 3);
			A.setZero();
			B.setZero();
			A.block(0, 0, 3, 3).setIdentity();
			Eigen::Vector3d st;
			st.setConstant(this->body_step_t);
			A.block(0, 3, 3, 3) = st.asDiagonal();
			A.block(3, 3, 3, 3).setIdentity();
			B.block(3, 0, 3, 3) = st.asDiagonal();
			traBodyGenerator.mpc_init(A, B, _Q, _F, _R, _W, 0);
		}
		void updateFootParams(const double& _t, const double& _balanceMarginX, const double& _balanceMarginY, const Eigen::Matrix<double, 16, 16>& _Q, const Eigen::Matrix<double, 16, 16>& _F, const Eigen::Matrix<double, 8, 8>& _R, const Eigen::Matrix<double, 8, 8>& _W)
		{
			this->foot_step_t = _t;
			this->balanceMarginX = _balanceMarginX;
			this->balanceMarginY = _balanceMarginY;
			MatrixXd A, B;
			A.resize(16, 16);
			B.resize(16, 8);
			A.setZero();
			B.setZero();
			A.block(0, 0, 8, 8).setIdentity();
			Eigen::Vector<double, 8> st;
			st.setConstant(this->foot_step_t);
			A.block(0, 8, 8, 8) = st.asDiagonal();
			A.block(8, 8, 8, 8).setIdentity();
			B.block(8, 0, 8, 8) = st.asDiagonal();
			traFootGenerator.mpc_init(A, B, _Q, _F, _R, _W, 0);
		}
		/* 更新等式约束 */
		void setBodyInputConstrain(const Vector3d& _aMax)
		{
			// 加速度上下限(对向量每个元素取绝对值操作)
			this->updateBodyIeqConstrain(_aMax);
			Eigen::Vector3d ub = _aMax.cwiseAbs();
			Eigen::Vector3d lb = -_aMax.cwiseAbs();
			traBodyGenerator.setConstrain(lb, ub);
		}
		void setFootInputConstrain(const Vector<double, 8>& _aMax)
		{
			Eigen::Vector<double, 8> ub = _aMax.cwiseAbs();
			Eigen::Vector<double, 8> lb = -_aMax.cwiseAbs();
			traFootGenerator.setConstrain(lb, ub);
		}
		/* 获取mpc原始输出 */
		const Eigen::Vector3d& getBodyMpcOriOut()
		{
			return this->bodyAcc;
		}
		const Eigen::Matrix<double, 3, 4>& getFootMpcOriOut()
		{
			return this->footAcc;
		}
	};

	class BalanceJointPlanning : public SupportPolygon
	{
	private:
		/* mpc轨迹规划器 */
		mpcCal<22, 11, 32, 1, 5> traGenerator;
		/* 每步预测的时间长度 */
		double step_t;
		/* mpc预测输出 */
		Eigen::Vector<double, 3> bodyAcc;
		Eigen::Matrix<double, 3, 4> footAcc;
		/* 坐标转换 */
		Eigen::Matrix3d Rsb_c;

		/* 配置关键不等式约束 */
		void updateJointIeqConstrain(const Eigen::Vector<double, 11>& _aMax)
		{
			// 平衡约束转化为线性锥约束
			Eigen::Matrix<double, 16, 3> cbA;
			cbA.setZero();
			cbA << 1, 0, this->balanceCone, -1, 0, this->balanceCone, 0, 1, this->balanceCone, 0, -1, this->balanceCone, \
				1, 1, 1.414 * this->balanceCone, -1, 1, 1.414 * this->balanceCone, -1, -1, 1.414 * this->balanceCone, 1, -1, 1.414 * this->balanceCone, \
				1, 1, -1.414 * this->balanceCone, -1, 1, -1.414 * this->balanceCone, -1, -1, -1.414 * this->balanceCone, 1, -1, 1.414 * this->balanceCone, \
				1, 1, 0, 1, -1, 0, -1, 1, 0, -1, -1, 0;
			Eigen::Vector<double, 16> Abub, Ablb;
			Abub.block(0, 0, 4, 1).setConstant(-(this->balanceCone * -9.81));
			Ablb.block(0, 0, 4, 1).setConstant(-100000.);
			Abub.block(4, 0, 4, 1).setConstant(-1.414 * this->balanceCone * -9.81);
			Ablb.block(4, 0, 4, 1).setConstant(-100000.);
			Abub.block(8, 0, 4, 1).setConstant(100000.);
			Ablb.block(8, 0, 4, 1).setConstant(1.414 * this->balanceCone * -9.81);
			Abub.block(12, 0, 4, 1).setConstant(_aMax.segment(0, 2).norm());
			Ablb.block(12, 0, 4, 1).setConstant(-_aMax.segment(0, 2).norm());

			// 对足端的加速度边界不等式约束设置所有边界
			Eigen::Matrix<double, 16, 11> cfA;
			cfA.setZero();
			Eigen::Vector<double, 16> Afub, Aflb;
			Afub.setConstant(100000.);
			for (int i = 0; i < 4; i++)
			{
				cfA(4*i + 0, 0) = -1; cfA(4*i + 0, 2) = -this->balanceMarginX;
				cfA(4*i + 1, 0) = 1;  cfA(4*i + 1, 2) = -this->balanceMarginX;
				cfA(4*i + 2, 1) = -1; cfA(4*i + 2, 2) = -this->balanceMarginY;
				cfA(4*i + 3, 1) = 1;  cfA(4*i + 3, 2) = -this->balanceMarginY;
				Aflb(4 * i + 0) = this->balanceMarginX * -9.81;
				Aflb(4 * i + 1) = this->balanceMarginX * -9.81;
				Aflb(4 * i + 2) = this->balanceMarginY * -9.81;
				Aflb(4 * i + 3) = this->balanceMarginY * -9.81;
			}
			cfA(0, 3) = 1; cfA(1, 3) = 1;
			cfA(2, 4) = 1; cfA(3, 4) = 1;
			cfA(4, 5) = 1; cfA(5, 5) = 1;
			cfA(6, 6) = -1; cfA(7, 6) = -1;
			cfA(8, 7) = -1; cfA(9, 7) = -1;
			cfA(10, 8) = 1; cfA(11, 8) = 1;
			cfA(12, 9) = -1; cfA(13, 9) = -1;
			cfA(14, 10) = -1; cfA(15, 10) = -1;

			// 联合约束矩阵
			Eigen::Matrix<double, 32, 11> cA;
			cA.setZero();
			cA.block(0, 0, 16, 3) = cbA;
			cA.block(16, 0, 16, 11) = cfA;
			Eigen::Vector<double, 32> Aub, Alb;
			Aub.segment(0, 16) = Abub;
			Alb.segment(0, 16) = Ablb;
			Aub.segment(16, 16) = Afub;
			Alb.segment(16, 16) = Aflb;

			traGenerator.setBoxConstrain(cA, Alb, Aub);
		}

		/* 求解规划器 */
		void getJointSolution()
		{
			static Eigen::Vector3d last_Vbr = this->Vbr;
			Eigen::Vector<double, 22> y, x;
			x.block(0, 0, 3, 1) = this->Pbr;
			x.block(3, 0, 3, 1) = this->Vbr;
			y.block(0, 0, 3, 1) = this->Pbt;
			y.block(3, 0, 3, 1) = this->Vbt;

			static Eigen::Matrix<double, 3, 4> last_Vlr = this->Vlr;
			x.block(6, 0, 2, 1) = this->Plr.col(0).segment(0, 2);
			x.block(8, 0, 2, 1) = this->Plr.col(1).segment(0, 2);
			x.block(10, 0, 2, 1) = this->Plr.col(2).segment(0, 2);
			x.block(12, 0, 2, 1) = this->Plr.col(3).segment(0, 2);
			x.block(14, 0, 2, 1) = this->Vlr.col(0).segment(0, 2);
			x.block(16, 0, 2, 1) = this->Vlr.col(1).segment(0, 2);
			x.block(18, 0, 2, 1) = this->Vlr.col(2).segment(0, 2);
			x.block(12, 0, 2, 1) = this->Vlr.col(3).segment(0, 2);
			y.block(6, 0, 2, 1) = this->Plt.col(0).segment(0, 2);
			y.block(8, 0, 2, 1) = this->Plt.col(1).segment(0, 2);
			y.block(10, 0, 2, 1) = this->Plt.col(2).segment(0, 2);
			y.block(12, 0, 2, 1) = this->Plt.col(3).segment(0, 2);
			y.block(14, 0, 8, 1).setZero();

			this->traGenerator.mpc_update(y, x, 100, 1.);
			this->traGenerator.mpc_solve();

			this->bodyAcc = this->traGenerator.getOutput().segment(0,3);
			this->Vbr = this->Vbr + 0.002 * this->Rsb_c * this->bodyAcc;
			this->Pbr = this->Pbr + 0.5 * 0.002 * (last_Vbr + this->Vbr);
			last_Vbr = this->Vbr;

			this->footAcc.block(0, 0, 2, 1) = this->traGenerator.getOutput().block(3, 0, 2, 1);
			this->footAcc.block(0, 1, 2, 1) = this->traGenerator.getOutput().block(5, 0, 2, 1);
			this->footAcc.block(0, 2, 2, 1) = this->traGenerator.getOutput().block(7, 0, 2, 1);
			this->footAcc.block(0, 3, 2, 1) = this->traGenerator.getOutput().block(9, 0, 2, 1);
			this->Vlr = this->Vlr + 0.002 * this->footAcc;
			this->Plr = this->Plr + 0.5 * 0.002 * (last_Vlr + this->Vlr);
			last_Vlr = this->Vlr;
		}

	public:
		/* 构造函数 */
		BalanceJointPlanning() : traGenerator(PL_LOW)
		{
			// 初始化参数
			this->Pbr.setZero();
			this->Vbr.setZero();
			this->Plr.setZero();
			this->Vlr.setZero();
			this->Rsb_c.setIdentity();
		}
		/* 规划器预热 */
		void warmUp(const Eigen::Matrix3d& _Rsbc)
		{
			this->Rsb_c = _Rsbc;
			this->updateSupportPolygon();
		}
		/* 使用机身规划器 */
		void useJointPlan(const Eigen::Vector<double,11>& _aMax)
		{
			this->updateJointIeqConstrain(_aMax);
			Eigen::Vector<double, 11> ub, lb;
			ub = _aMax.cwiseAbs();
			lb = -_aMax.cwiseAbs();
			this->traGenerator.setConstrain(lb, ub);
			this->getJointSolution();
		}
		void updateJointParams(double _t, double _bMargin, double _bMarginX, double _bMarginY, const Eigen::Matrix<double, 22, 22>& _Q, const Eigen::Matrix<double, 22, 22>& _F, const Eigen::Matrix<double, 11, 11>& _R, const Eigen::Matrix<double, 11, 11>& _W)
		{
			this->step_t = _t;
			this->balanceMargin = _bMargin;
			this->balanceMarginX = _bMarginX;
			this->balanceMarginY = _bMarginY;
			MatrixXd A, B;
			A.resize(22, 22);
			B.resize(22, 11);
			A.setZero();
			B.setZero();
			Eigen::Vector<double, 3> bt;
			bt.setConstant(this->step_t);
			A.block(0, 0, 3, 3).setIdentity();
			A.block(0, 3, 3, 3) = bt.asDiagonal();
			A.block(3, 3, 3, 3).setIdentity();
			B.block(3, 0, 3, 3) = bt.asDiagonal();
			Eigen::Vector<double, 8> ft;
			ft.setConstant(this->step_t);
			A.block(6, 6, 8, 8).setIdentity();
			A.block(6, 14, 8, 8) = ft.asDiagonal();
			A.block(14, 14, 8, 8).setIdentity();
			B.block(14, 3, 8, 8) = ft.asDiagonal();
			traGenerator.mpc_init(A, B, _Q, _F, _R, _W, 0);
		}
	};

}