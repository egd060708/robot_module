#pragma once

#include <Eigen/Dense>
#include "mpcCal.h"
#include <iostream>
#include "mathTool.h"

using namespace Eigen;

namespace Quadruped {

	class PlanBase {
	protected:
		/* 机身的轨迹规划 */
		Vector3d Pbt,Pbc,Vbt,Vbc,Pwc,Vwc;
		/* 足端的轨迹规划 */
		Eigen::Matrix<double, 3, 4> Plt, Plc, Vlt, Vlc;
		/* 规划结果 */
		Vector3d Pbr, Vbr;
		Eigen::Matrix<double, 3, 4> Plr, Vlr, Pli;
		/* 坐标转换 */
		Eigen::Matrix3d Rsb_c;
		Eigen::Matrix3d Rsbh_c;
		Eigen::Matrix4d Tsb_c;
		Eigen::Matrix4d Tsbh_c;
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
		void updateRsb(const Matrix3d& _Rsb, const Matrix3d& _Rsbh, const Matrix4d& _Tsb, const Matrix4d& _Tsbh)
		{
			this->Rsb_c = _Rsb;
			this->Rsbh_c = _Rsbh;
			this->Tsb_c = _Tsb;
			this->Tsbh_c = _Tsbh;
		}
		void updateBodyHighLevelTar(const Vector3d& _Pbt, const Vector3d& _Vbt)
		{
			this->Pbt = _Pbt;
			this->Vbt = _Vbt;
			//std::cout << "tar: " << _Pbt << std::endl;
		}
		void updateFootHighLevelTar(const Eigen::Matrix<double, 3, 4>& _Plt, const Eigen::Matrix<double, 3, 4>& _Vlt)
		{
			/*for (int i = 0; i < 4; i++)
			{
				this->Plt.col(i) = this->Rsb_c.transpose() * _Plt.col(i);
			}*/
			this->Plt = _Plt;
			this->Vlt = _Vlt;
		}
		void updateBodyState(const Vector3d& _Pbc, const Vector3d& _Vbc)
		{
			this->Pbc = _Pbc;
			this->Vbc = _Vbc;
			//std::cout << "cur: " << _Pbc << std::endl;
		}
		void updateFootState(const Eigen::Matrix<double, 3, 4>& _Plc, const Eigen::Matrix<double, 3, 4>& _Vlc)
		{
			this->Plc = _Plc;
			this->Vlc = _Vlc;
		}
		void updateWBodyState(const Vector3d& _Pwc, const Vector3d& _Vwc)
		{
			this->Pwc = _Pwc;
			this->Vwc = _Vwc;
		}
		virtual Vector3d getBodyPlanPosition()
		{
			return this->Pbr;
		}
		virtual Vector3d getBodyPlanVelocity()
		{
			return this->Vbr;
		}
		virtual Eigen::Matrix<double, 3, 4> getFootPlanPosition()
		{
			return this->Plr;
		}
		virtual Eigen::Matrix<double, 3, 4> getFootPlanVelocity()
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
			this->Pli = _Plr;
			this->Plr = _Plr;
		}
		void setInitFootPlanVelocity(const Eigen::Matrix<double, 3, 4>& _Vlr)
		{
			this->Vlr = _Vlr;
		}

		Eigen::VectorXd normConstrain(const Eigen::VectorXd& _cur, const Eigen::VectorXd& _ref, const double& _value)
		{
			double scale = (_cur - _ref).norm() / _value;
			Eigen::VectorXd ret;
			if (scale > 1.)
			{
				ret = _ref + (_cur - _ref).normalized() * _value;
			}
			else
			{
				ret = _cur;
			}
			return ret;
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
			double temp = -(this->Plc.col(1) - this->Plc.col(0)).transpose() * (this->Plc.col(0) - this->Pbc);
			tumblingAxis(0) = temp / (this->Plc.col(1) - this->Plc.col(0)).squaredNorm();
			tumblingPoint.col(0) = this->Plc.col(0) + tumblingAxis(0) * (this->Plc.col(1) - this->Plc.col(0));

			temp = -(this->Plc.col(3) - this->Plc.col(1)).transpose() * (this->Plc.col(1) - this->Pbc);
			tumblingAxis(1) = temp / (this->Plc.col(3) - this->Plc.col(1)).squaredNorm();
			tumblingPoint.col(1) = this->Plc.col(1) + tumblingAxis(1) * (this->Plc.col(3) - this->Plc.col(1));

			temp = -(this->Plc.col(2) - this->Plc.col(3)).transpose() * (this->Plc.col(3) - this->Pbc);
			tumblingAxis(2) = temp / (this->Plc.col(2) - this->Plc.col(3)).squaredNorm();
			tumblingPoint.col(2) = this->Plc.col(3) + tumblingAxis(2) * (this->Plc.col(2) - this->Plc.col(3));

			temp = -(this->Plc.col(0) - this->Plc.col(2)).transpose() * (this->Plc.col(2) - this->Pbc);
			tumblingAxis(3) = temp / (this->Plc.col(0) - this->Plc.col(2)).squaredNorm();
			tumblingPoint.col(3) = this->Plc.col(2) + tumblingAxis(3) * (this->Plc.col(0) - this->Plc.col(2));
		}

		/* 找到平衡锥约束参数 */
		void findBalanceCone()
		{
			Eigen::Vector3d g(0, 0, -9.81);
			// 找到最严格的约束点(计算最小夹角)，以及相对于x轴和y轴的最严格约束点
			double temp = g.dot(this->tumblingPoint.col(0) - this->Pbc) / (g.norm() * (this->tumblingPoint.col(0) - this->Pbc).norm());
			double tempx = temp;
			double tempy = g.dot(this->tumblingPoint.col(1) - this->Pbc) / (g.norm() * (this->tumblingPoint.col(1) - this->Pbc).norm());
			int max_idx = 0, max_idxX = 0, max_idxY = 1;
			Eigen::Vector4d temp4;
			for (int i = 0; i < 4; i++)
			{
				double temp1 = g.dot(this->tumblingPoint.col(i) - this->Pbc) / (g.norm() * (this->tumblingPoint.col(i) - this->Pbc).norm());
				if (temp1 > temp)
				{
					max_idx = i;
					temp = temp1;
				}
				if (i == 0 || i == 2)
				{
					if (temp1 > tempx)
					{
						max_idxX = i;
						tempx = temp1;
					}
				}
				else if (i == 1 || i == 3)
				{
					if (temp1 > tempy)
					{
						max_idxY = i;
						tempy = temp1;
					}
				}
			}
			balanceCone = balanceMargin * (this->tumblingPoint.col(max_idx) - this->Pbc).block(0, 0, 2, 1).norm() / abs((this->tumblingPoint.col(max_idx) - this->Pbc)(2));
			balanceConeX = balanceMarginX * (this->tumblingPoint.col(max_idxX) - this->Pbc).block(0, 0, 2, 1).norm() / abs((this->tumblingPoint.col(max_idxX) - this->Pbc)(2));
			balanceConeY = balanceMarginY * (this->tumblingPoint.col(max_idxY) - this->Pbc).block(0, 0, 2, 1).norm() / abs((this->tumblingPoint.col(max_idxY) - this->Pbc)(2));
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
		mpcCal<6, 3, 8, 1, 10> traBodyGenerator;
		/* 预测时间步长 */
		double body_step_t = 0;
		double foot_step_t = 0;
		/* 机身mpc输出观测 */
		Eigen::Vector3d bodyAcc;
		/* 针对足端的规划生成器 */
		mpcCal<16, 8, 16, 1, 10> traFootGenerator;
		/* 足端mpc输出观测 */
		Eigen::Matrix<double, 3, 4> footAcc;

		/* 更新机身不等式约束 */
		void updateBodyIeqConstrain(const Vector3d& _aMax)
		{
			// 平衡约束转化为线性锥约束
			Eigen::Matrix<double, 8, 3> cA;
			cA.setZero();
			cA << 1, 0, this->balanceConeX, -1, 0, this->balanceConeX, 0, 1, this->balanceConeY, 0, -1, this->balanceConeY, \
				/*1, 1, 1.414 * this->balanceCone, -1, 1, 1.414 * this->balanceCone, -1, -1, 1.414 * this->balanceCone, 1, -1, 1.414 * this->balanceCone, \
				1, 1, -1.414 * this->balanceCone, -1, 1, -1.414 * this->balanceCone, -1, -1, -1.414 * this->balanceCone, 1, -1, 1.414 * this->balanceCone, \*/
				1, 1, 0, 1, -1, 0, -1, 1, 0, -1, -1, 0;
			Eigen::Vector<double, 8> Aub, Alb;
			Aub.block(0, 0, 2, 1).setConstant(-(this->balanceConeX * -9.81));
			Aub.block(2, 0, 2, 1).setConstant(-(this->balanceConeY * -9.81));
			Alb.block(0, 0, 4, 1).setConstant(-100000.);
			//Aub.block(4, 0, 4, 1).setConstant(-1.414 * this->balanceCone * -9.81);
			//Alb.block(4, 0, 4, 1).setConstant(-100000.);
			//Aub.block(8, 0, 4, 1).setConstant(100000.);
			//Alb.block(8, 0, 4, 1).setConstant(1.414 * this->balanceCone * -9.81);
			Aub.block(4, 0, 4, 1).setConstant(_aMax.segment(0, 2).norm());
			Alb.block(4, 0, 4, 1).setConstant(-_aMax.segment(0, 2).norm());
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
			this->traBodyGenerator.mpc_solve(1);
			this->bodyAcc = this->traBodyGenerator.getOutput();
			this->Vbr = this->Vbr + 0.002 * this->Rsbh_c * this->bodyAcc;
			this->Pbr = this->Pbr + 0.5 * 0.002 * (last_Vbr + this->Vbr);
			last_Vbr = this->Vbr;
		}

		/* 更新足端不等式约束 */
		void updateFootIeqConstrain()
		{
			// 对足端的加速度边界不等式约束设置所有边界
			double edgeX1, edgeX2, edgeY1, edgeY2;
			Eigen::Vector3d realBodyAcc = this->bodyAcc;
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
			this->traFootGenerator.mpc_solve(1);
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
			this->Rsbh_c.setIdentity();
		}
		/* 规划器预热 */
		void warmUp()
		{
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
			B.block(3, 0, 3, 3) = st.asDiagonal()*this->Rsbh_c;
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
		mpcCal<22, 11, 24, 1, 8> traGenerator;
		/* 每步预测的时间长度 */
		double step_t;
		/* mpc预测输出 */
		Eigen::Vector<double, 3> bodyAcc;
		Eigen::Matrix<double, 3, 4> footAcc;
		/* 足端加速度约束与机身稳定裕度的比例 */
		double fbp;

		/* 配置关键不等式约束 */
		void updateJointIeqConstrain(const Eigen::Vector<double, 11>& _aMax, const Eigen::Vector3d& _bhAcc)
		{
			// 平衡约束转化为线性锥约束
			Eigen::Matrix<double, 8, 3> cbA;
			cbA.setZero();
			cbA << 1, 0, this->balanceConeX, -1, 0, this->balanceConeX, 0, 1, this->balanceConeY, 0, -1, this->balanceConeY, \
				1, 1, 0, 1, -1, 0, -1, 1, 0, -1, -1, 0;
			/*cbA <<0, 0, this->balanceConeX, 0, 0, this->balanceConeY, 0, 0, this->balanceConeX, 0, 0, this->balanceConeY, \
				1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1, 0;*/
			Eigen::Vector<double, 8> Abub, Ablb;
			Abub.block(0, 0, 2, 1).setConstant(-(this->balanceConeX * -9.81));
			Abub.block(2, 0, 2, 1).setConstant(-(this->balanceConeY * -9.81));
			Abub.block(4, 0, 4, 1).setConstant(_aMax.segment(0, 2).norm());
			/*Abub.block(0, 0, 2, 1).setConstant(-abs(_bhAcc(0)) - (this->balanceConeX * -9.81));
			Abub.block(2, 0, 2, 1).setConstant(-abs(_bhAcc(1)) - (this->balanceConeY * -9.81));
			Abub.block(4, 0, 2, 1).setConstant(-this->balanceConeX * (_bhAcc(2) - 2 * 9.81));
			Abub.block(6, 0, 2, 1).setConstant(-this->balanceConeY * (_bhAcc(2) - 2 * 9.81));*/

			Ablb.block(0, 0, 4, 1).setConstant(-100000.);
			Ablb.block(4, 0, 4, 1).setConstant(-_aMax.segment(0, 2).norm());
			//Ablb.block(0, 0, 8, 1).setConstant(-100000.);


			// 对足端的加速度边界不等式约束设置所有边界
			Eigen::Matrix<double, 16, 11> cfA;
			cfA.setZero();
			Eigen::Vector<double, 16> Afub, Aflb;
			Afub.setConstant(100000.);
			for (int i = 0; i < 4; i++)
			{
				// 机身x和z加速度受影响
				cfA(4*i + 0, 0) = -1; cfA(4*i + 0, 2) = -this->fbp * this->balanceConeX;
				cfA(4*i + 1, 0) = 1;  cfA(4*i + 1, 2) = -this->fbp * this->balanceConeX;
				cfA(4*i + 2, 1) = -1; cfA(4*i + 2, 2) = -this->fbp * this->balanceConeY;
				cfA(4*i + 3, 1) = 1;  cfA(4*i + 3, 2) = -this->fbp * this->balanceConeY;
				Aflb(4 * i + 0) = this->fbp * this->balanceConeX * (-9.81);
				Aflb(4 * i + 1) = this->fbp * this->balanceConeX * (-9.81);
				Aflb(4 * i + 2) = this->fbp * this->balanceConeY * (-9.81);
				Aflb(4 * i + 3) = this->fbp * this->balanceConeY * (-9.81);
				// 机身水平方向加速度受影响
				/*cfA(4 * i + 0, 0) = -1;
				cfA(4 * i + 1, 0) = 1; 
				cfA(4 * i + 2, 1) = -1;
				cfA(4 * i + 3, 1) = 1; 
				Aflb(4 * i + 0) = this->fbp * this->balanceConeX * (this->bodyAcc(2)-9.81);
				Aflb(4 * i + 1) = this->fbp * this->balanceConeX * (this->bodyAcc(2)-9.81);
				Aflb(4 * i + 2) = this->fbp * this->balanceConeY * (this->bodyAcc(2)-9.81);
				Aflb(4 * i + 3) = this->fbp * this->balanceConeY * (this->bodyAcc(2)-9.81);*/
				// 机身加速度不受足端影响
				/*Aflb(4 * i + 0) = this->bodyAcc(0) + this->fbp * this->balanceConeX * (this->bodyAcc(2) - 9.81);
				Aflb(4 * i + 1) = -this->bodyAcc(0) + this->fbp * this->balanceConeX * (this->bodyAcc(2) - 9.81);
				Aflb(4 * i + 2) = this->bodyAcc(1) + this->fbp * this->balanceConeY * (this->bodyAcc(2) - 9.81);
				Aflb(4 * i + 3) = -this->bodyAcc(1) + this->fbp * this->balanceConeY * (this->bodyAcc(2) - 9.81);*/
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
			Eigen::Matrix<double, 24, 11> cA;
			cA.setZero();
			cA.block(0, 0, 8, 3) = cbA;
			cA.block(8, 0, 16, 11) = cfA;
			Eigen::Vector<double, 24> Aub, Alb;
			Aub.segment(0, 8) = Abub;
			Alb.segment(0, 8) = Ablb;
			Aub.segment(8, 16) = Afub;
			Alb.segment(8, 16) = Aflb;

			/*std::cout << "cA: \n" << cA << std::endl;
			std::cout << "Aub:\n" << Aub.transpose() << std::endl;
			std::cout << "Alb:\n" << Alb.transpose() << std::endl;*/

			traGenerator.setBoxConstrain(cA, Alb, Aub);

			// 根据加速度观测增加额外的约束
			static MeanFilter<10> accF[3];
			double bhAccX = accF[0].f(_bhAcc(0));
			double bhAccY = accF[1].f(_bhAcc(1));
			double bhAccZ = accF[2].f(_bhAcc(2));
			Eigen::Vector<double, 11> ub, lb;
			ub = _aMax.cwiseAbs();
			lb = -_aMax.cwiseAbs();
			if ((-abs(bhAccX) - (this->balanceConeX * -9.81)) / this->balanceConeX < (-abs(bhAccY) - (this->balanceConeY * -9.81)) / this->balanceConeY)
			{
				ub(2) = (-abs(bhAccX) - (this->balanceConeX * -9.81)) / this->balanceConeX;
			}
			else
			{
				ub(2) = (-abs(bhAccY) - (this->balanceConeY * -9.81)) / this->balanceConeY;
			}
			ub(2) = constrain(ub(2), _aMax.cwiseAbs()(2), -_aMax.cwiseAbs()(2));
			if (bhAccZ > 18.)
			{
				bhAccZ = 18.;
			}
			ub(0) = constrain(-this->balanceConeX * (bhAccZ - 2 * 9.81), _aMax.cwiseAbs()(0), -_aMax.cwiseAbs()(0));
			lb(0) = constrain(this->balanceConeX * (bhAccZ - 2 * 9.81), _aMax.cwiseAbs()(0), -_aMax.cwiseAbs()(0));
			ub(1) = constrain(-this->balanceConeY * (bhAccZ - 2 * 9.81), _aMax.cwiseAbs()(1), -_aMax.cwiseAbs()(1));
			lb(1) = constrain(this->balanceConeY * (bhAccZ - 2 * 9.81), _aMax.cwiseAbs()(1), -_aMax.cwiseAbs()(1));

			double edge[4] = { -bhAccX + this->fbp * this->balanceConeX * (bhAccZ - 2 * 9.81), \
								bhAccX + this->fbp * this->balanceConeX * (bhAccZ - 2 * 9.81), \
								-bhAccY + this->fbp * this->balanceConeX * (bhAccZ - 2 * 9.81), \
								bhAccY + this->fbp * this->balanceConeX * (bhAccZ - 2 * 9.81) };
			lb(3) = bhAccX > 0 ? edge[0] : edge[1];
			lb(4) = bhAccY > 0 ? edge[2] : edge[3];
			lb(5) = bhAccX > 0 ? edge[0] : edge[1];
			ub(6) = -(bhAccY > 0 ? edge[2] : edge[3]);
			ub(7) = -(bhAccX > 0 ? edge[0] : edge[1]);
			lb(8) = bhAccY > 0 ? edge[2] : edge[3];
			ub(9) = -(bhAccX > 0 ? edge[0] : edge[1]);
			ub(10) = -(bhAccY > 0 ? edge[2] : edge[3]);

			lb(3) = constrain(lb(3), _aMax.cwiseAbs()(3), -_aMax.cwiseAbs()(3));
			lb(4) = constrain(lb(4), _aMax.cwiseAbs()(4), -_aMax.cwiseAbs()(4));
			lb(5) = constrain(lb(5), _aMax.cwiseAbs()(5), -_aMax.cwiseAbs()(5));
			ub(6) = constrain(ub(6), _aMax.cwiseAbs()(6), -_aMax.cwiseAbs()(6));
			ub(7) = constrain(ub(7), _aMax.cwiseAbs()(7), -_aMax.cwiseAbs()(7));
			lb(8) = constrain(lb(8), _aMax.cwiseAbs()(8), -_aMax.cwiseAbs()(8));
			ub(9) = constrain(ub(9), _aMax.cwiseAbs()(9), -_aMax.cwiseAbs()(9));
			ub(10) = constrain(ub(10), _aMax.cwiseAbs()(10), -_aMax.cwiseAbs()(10));
			this->traGenerator.setConstrain(lb, ub);
			//std::cout << "ub: " << ub.block(0, 0, 3, 1).transpose() << std::endl;
			//std::cout << "lb: " << lb.block(0, 0, 3, 1).transpose() << std::endl;
		}

		/* 求解规划器 */
		void getJointSolution()
		{
			static Eigen::Vector3d last_Vbr = this->Vbr;
			Eigen::Vector<double, 22> y, x;
			Eigen::Vector3d Py(this->Pbr(0), this->Pbr(1), this->Pbt(2));
			x.block(0, 0, 3, 1) = this->Pbr;
			x.block(3, 0, 3, 1) = this->Vbr;
			//y.block(0, 0, 3, 1) = this->Pbt;
			y.block(0, 0, 3, 1) = Py;
			y.block(3, 0, 3, 1) = this->Vbt;

			static Eigen::Matrix<double, 3, 4> last_Vlr = this->Vlr;
			x.block(6, 0, 2, 1) = this->Plr.col(0).segment(0, 2);
			x.block(8, 0, 2, 1) = this->Plr.col(1).segment(0, 2);
			x.block(10, 0, 2, 1) = this->Plr.col(2).segment(0, 2);
			x.block(12, 0, 2, 1) = this->Plr.col(3).segment(0, 2);
			x.block(14, 0, 2, 1) = this->Vlr.col(0).segment(0, 2);
			x.block(16, 0, 2, 1) = this->Vlr.col(1).segment(0, 2);
			x.block(18, 0, 2, 1) = this->Vlr.col(2).segment(0, 2);
			x.block(20, 0, 2, 1) = this->Vlr.col(3).segment(0, 2);
			y.block(6, 0, 2, 1) = this->Plt.col(0).segment(0, 2);
			y.block(8, 0, 2, 1) = this->Plt.col(1).segment(0, 2);
			y.block(10, 0, 2, 1) = this->Plt.col(2).segment(0, 2);
			y.block(12, 0, 2, 1) = this->Plt.col(3).segment(0, 2);
			y.block(14, 0, 8, 1).setZero();

			this->traGenerator.mpc_update(y, x, 1000, 0.02);
			this->traGenerator.mpc_solve(0);

			this->bodyAcc = this->traGenerator.getOutput().block(0,0,3,1);
			// 生成机身速度
			this->Vbr = this->Vbr + 0.002 * this->Rsbh_c * this->bodyAcc;
			// 速度约束
			/*this->Vbr(0) = slopeConstrain(this->Vbr(0), this->Vbc(0), 0.1, -0.1);
			this->Vbr(1) = slopeConstrain(this->Vbr(1), this->Vbc(1), 0.1, -0.1);*/
			/*this->Vbr(0) = constrain(this->Vbr(0), 2.5, -2.5);
			this->Vbr(1) = constrain(this->Vbr(1), 2.5, -2.5);*/
			this->Vbr.block(0, 0, 2, 1) = this->normConstrain(this->Vbr.block(0, 0, 2, 1), this->Vbc.block(0, 0, 2, 1), 1.5);
			this->Vbr(2) = constrain(this->Vbr(2), 0.5, -0.5);
			// 生成机身位置
			this->Pbr = this->Pbr + 0.5 * 0.002 * (last_Vbr + this->Vbr);
			// 位置约束
			/*this->Pbr(0) = slopeConstrain(this->Pbr(0), this->Pbc(0), 0.3, -0.3);
			this->Pbr(1) = slopeConstrain(this->Pbr(1), this->Pbc(1), 0.3, -0.3);*/
			//this->Pbr(2) = slopeConstrain(this->Pbr(2), this->Pbc(2), 0.01, -0.01);
			this->Pbr.block(0, 0, 2, 1) = this->normConstrain(this->Pbr.block(0, 0, 2, 1), this->Pbc.block(0, 0, 2, 1), 0.3);
			this->Pbr(2) = constrain(this->Pbr(2), 0.6, 0.45);
			last_Vbr = this->Vbr;

			//std::cout << "gen: " << this->Pbr << std::endl;

			this->footAcc.block(0, 0, 2, 1) = this->traGenerator.getOutput().block(3, 0, 2, 1);
			this->footAcc.block(0, 1, 2, 1) = this->traGenerator.getOutput().block(5, 0, 2, 1);
			this->footAcc.block(0, 2, 2, 1) = this->traGenerator.getOutput().block(7, 0, 2, 1);
			this->footAcc.block(0, 3, 2, 1) = this->traGenerator.getOutput().block(9, 0, 2, 1);
			// 生成接触点速度
			this->Vlr = this->Vlr + 0.002 * this->footAcc;
			// 生成接触点位置
			this->Plr = this->Plr + 0.5 * 0.002 * (last_Vlr + this->Vlr);
			last_Vlr = this->Vlr;
			// 位置约束
			for (int i = 0; i < 4; i++)
			{
				this->Plr.col(i) = this->normConstrain(this->Plr.col(i), this->Pli.col(i), 0.8);
			}

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
			this->Rsbh_c.setIdentity();
			this->Tsb_c.setZero();
			this->Tsb_c.block(0, 0, 3, 3).setIdentity();
			this->Tsb_c(3, 3) = 1;
			this->Tsbh_c.setZero();
			this->Tsbh_c.block(0, 0, 3, 3).setIdentity();
			this->Tsbh_c(3, 3) = 1;
			fbp = 1.;
		}
		/* 规划器预热 */
		void warmUp()
		{
			this->updateSupportPolygon();
		}
		/* 使用机身规划器 */
		void useJointPlan(const Eigen::Vector<double,11>& _aMax, const Eigen::Vector3d& _bhAcc)
		{
			this->updateJointIeqConstrain(_aMax, _bhAcc);
			this->getJointSolution();
		}
		void updateJointParams(double _t, double _fbp, double _bMarginX, double _bMarginY, const Eigen::Matrix<double, 22, 22>& _Q, const Eigen::Matrix<double, 22, 22>& _F, const Eigen::Matrix<double, 11, 11>& _R, const Eigen::Matrix<double, 11, 11>& _W)
		{
			this->step_t = _t;
			this->fbp = _fbp;
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
			B.block(3, 0, 3, 3) = bt.asDiagonal() * this->Rsbh_c;
			Eigen::Vector<double, 8> ft;
			ft.setConstant(this->step_t);
			A.block(6, 6, 8, 8).setIdentity();
			A.block(6, 14, 8, 8) = ft.asDiagonal();
			A.block(14, 14, 8, 8).setIdentity();
			B.block(14, 3, 8, 8) = ft.asDiagonal();
			traGenerator.mpc_init(A, B, _Q, _F, _R, _W, 0);
		}
		Eigen::Vector3d getBodyPlanPosition() override
		{
			return this->Pbr;
		}
		Eigen::Vector3d getBodyPlanVelocity() override
		{
			return this->Vbr;
		}
		Eigen::Matrix<double, 3, 4> getFootPlanPosition() override
		{
			return this->Plr;
		}
		Eigen::Matrix<double, 3, 4> getFootPlanVelocity() override
		{
			Eigen::Matrix<double, 3, 4> fVbr;
			fVbr.col(0) = this->Rsbh_c.transpose() * this->Vbr;
			fVbr.col(1) = this->Rsbh_c.transpose() * this->Vbr;
			fVbr.col(2) = this->Rsbh_c.transpose() * this->Vbr;
			fVbr.col(3) = this->Rsbh_c.transpose() * this->Vbr;
			//std::cout << fVbr << std::endl;
			return this->Vlr + fVbr;
		}
		Eigen::Vector3d getBodyPlanVelocityB()
		{
			return this->Rsb_c * this->Vbr;
		}
		Eigen::Matrix<double, 3, 4> getFootPlanPositionWorld()
		{
			Eigen::Matrix<double, 3, 4> realPlr;
			Eigen::Vector4d Pbi(0, 0, 0, 1);
			Pbi.block(0, 0, 3, 1) = this->Pbr;
			Pbi = this->Tsbh_c.inverse() * Pbi;
			for (int i = 0; i < 4; i++)
			{
				realPlr.col(i) = Pbi.block(0,0,3,1) + this->Plr.col(i);
				realPlr.col(i)(2) = 0;
			}
			return realPlr;
		}
	};

}