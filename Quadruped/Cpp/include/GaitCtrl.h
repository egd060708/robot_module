#pragma once
#include <math.h>
#include <iostream>
#include "BodyCtrlNew.h"
#include "LegCtrl.h"
#include "mathTool.h"
#include <vector>
using namespace std;

#define M_PI 3.14159265358979323846

namespace Quadruped
{
	enum class WaveStatus {
		STANCE_ALL = 0,
		SWING_ALL = 1,
		WAVE_ALL = 2,
		ADAPT = 3
	};

	class GaitCtrl
	{
	public:
		GaitCtrl(CtrlBase* _bc, LegCtrl* _lc[4], int timeStep, Eigen::Vector4d* _gaitPhase, Eigen::Vector4i* _gaitContact);
		Eigen::Vector3d calFootPos(int legID, Eigen::Vector2d vxyTargetGlobal, double dYawTarget, double phase);
		Eigen::Vector3d calFootPosW(int legID, Eigen::Vector2d vxyTargetGlobal, double dYawTarget, double phase, double maxStepL,const Eigen::Matrix3d& _slope);
		void initExpectK(Eigen::Vector3d _k);// 初始化期望运动控制参数
		void _updateFootPoints();// 更新中性立足点
		void initSwingParams(double _period, double _stancePhaseRatio, Eigen::Vector4d _bias, double _t);// 初始化摆动相关参数
		void calcWave(Eigen::Vector4d& phase, Eigen::Vector4i& contact, WaveStatus status, double _t, Eigen::Vector4d& _estPhase, Eigen::Vector4i& _estContact);
		void calcContactPhase(WaveStatus status, double _t, Eigen::Vector4d& _estPhase, Eigen::Vector4i& _estContact);

		void footPointMarking();
		void footStateTransform(double _t);
		void footAdapt();
		
		void setGait(Eigen::Vector2d vxyTargetGlobal, double dYawTarget, double gaitHeight);
		void run(Eigen::Matrix<double, 3, 4>& _feetPos, Eigen::Matrix<double, 3, 4>& _feetVel, double _maxStepL, const Eigen::Matrix3d& _slope);
		Eigen::Vector3d getFootPos(int i);
		Eigen::Vector3d getFootVel(int i);
		void restart();
	public:
		double dt;
		// 联合的控制器,机身控制器以及腿部控制器
		CtrlBase* bodyController;
		LegCtrl* legController[4];
		// 末端落脚点计算部分
		Eigen::Vector3d nextStep = Eigen::Vector3d::Zero();
		Eigen::Vector3d footPos = Eigen::Vector3d::Zero();
		Eigen::Vector3d bodyVelGlobal = Eigen::Vector3d::Zero();// 在世界坐标系下的机体线速度
		Eigen::Vector3d bodyAccGlobal = Eigen::Vector3d::Zero();// 在世界坐标系下的机体线加速度
		Eigen::Vector3d bodyWGlobal = Eigen::Vector3d::Zero();// 在世界坐标系下的机体角速度
		Eigen::Vector4d feetRadius = Eigen::Vector4d::Zero();// 足端运动半径
		Eigen::Vector4d feetInitAngle = Eigen::Vector4d::Zero();;// 足端初始运动角度
		double yaw, dYaw, nextYaw;
		double Tstance, Tswing;
		double kx, ky, kyaw;
		Eigen::Vector3d initLeg[4];
		Eigen::Vector4d fpmark = Eigen::Vector4d::Ones();
		Eigen::Vector4i recoverflag = Eigen::Vector4i::Zero();
		Eigen::Vector4i recoverflaglast = Eigen::Vector4i::Zero();
		// 接触状态与摆动状态检测
		double period;// 步态周期p
		double stRatio;// 触地系数r(归一化)
		Eigen::Vector4d estStRatio;// 用于被动预测的触地系数
		Eigen::Vector4d bias = Eigen::Vector4d::Zero();// 步态偏移系数，单腿偏移时间与步态周期的比值(归一化)
		Eigen::Vector4d normalT = Eigen::Vector4d::Zero();// 归一化的时间
		Eigen::Vector4d phase = Eigen::Vector4d::Zero();
		Eigen::Vector4d phasePast = Eigen::Vector4d::Zero();
		Eigen::Vector4i contact = Eigen::Vector4i::Zero();
		Eigen::Vector4i contactPast = Eigen::Vector4i::Zero();
		Eigen::Vector4i switchStatus = Eigen::Vector4i::Zero();// 是否能够切换状态，1为可以，0不可
		WaveStatus statusPast;
		Eigen::Vector4d passT;
		Eigen::Vector4d startT;
		// 摆动曲线生成
		double gaitHeight;
		Eigen::Vector2d VxyTarget;
		double dYawTarget;
		Eigen::Vector4d* gaitPhase;
		Eigen::Vector4d gaitPhasePast;
		Eigen::Vector4i* gaitContact;
		Eigen::Matrix<double, 3, 4> startP, endP, idealP, pastP;
		bool is_firstRun;
		// 腿部摆线生成函数
		double cycloidXYPosition(double startXY, double endXY, double phase);
		double cycloidXYVelocity(double startXY, double endXY, double phase);
		double cycloidZPosition(double startZ, double height, double phase);
		double cycloidZVelocity(double height, double phase);
	};

	GaitCtrl::GaitCtrl(CtrlBase* _bc, LegCtrl* _lc[4], int timeStep, Eigen::Vector4d* _gaitPhase, Eigen::Vector4i* _gaitContact)
	{
		this->gaitPhase = _gaitPhase;
		this->gaitContact = _gaitContact;
		this->bodyController = _bc;
		for (int i = 0; i < 4; i++)
		{
			this->legController[i] = _lc[i];
			initLeg[i].setZero();
		}
		
		initLeg[RB](0) = bodyController->bodyObject->initLegsXYPosition(0,RB);
		initLeg[RB](1) = bodyController->bodyObject->initLegsXYPosition(1,RB);
		initLeg[LB](0) = bodyController->bodyObject->initLegsXYPosition(0,LB);
		initLeg[LB](1) = bodyController->bodyObject->initLegsXYPosition(1,LB);
		initLeg[RF](0) = bodyController->bodyObject->initLegsXYPosition(0,RF);
		initLeg[RF](1) = bodyController->bodyObject->initLegsXYPosition(1,RF);
		initLeg[LF](0) = bodyController->bodyObject->initLegsXYPosition(0,LF);
		initLeg[LF](1) = bodyController->bodyObject->initLegsXYPosition(1,LF);
		for (int i = 0; i < 4; i++)
		{
			feetRadius(i) = sqrt(pow(initLeg[i](0), 2) + pow(initLeg[i](1), 2));
			feetInitAngle(i) = atan2(initLeg[i](1), initLeg[i](0));
		}
		dt = 0.001 * timeStep;
	}

	void GaitCtrl::initExpectK(Eigen::Vector3d _k)
	{
		kx = _k(0);
		ky = _k(1);
		kyaw = _k(2);
	}

	void GaitCtrl::_updateFootPoints()
	{
		initLeg[RB](0) = bodyController->bodyObject->initLegsXYPosition(0, RB);
		initLeg[RB](1) = bodyController->bodyObject->initLegsXYPosition(1, RB);
		initLeg[LB](0) = bodyController->bodyObject->initLegsXYPosition(0, LB);
		initLeg[LB](1) = bodyController->bodyObject->initLegsXYPosition(1, LB);
		initLeg[RF](0) = bodyController->bodyObject->initLegsXYPosition(0, RF);
		initLeg[RF](1) = bodyController->bodyObject->initLegsXYPosition(1, RF);
		initLeg[LF](0) = bodyController->bodyObject->initLegsXYPosition(0, LF);
		initLeg[LF](1) = bodyController->bodyObject->initLegsXYPosition(1, LF);
		for (int i = 0; i < 4; i++)
		{
			feetRadius(i) = sqrt(pow(initLeg[i](0), 2) + pow(initLeg[i](1), 2));
			feetInitAngle(i) = atan2(initLeg[i](1), initLeg[i](0));
		}
	}

	void GaitCtrl::initSwingParams(double _period, double _stancePhaseRatio, Eigen::Vector4d _bias, double _t)
	{
		period = _period;
		stRatio = _stancePhaseRatio;
		estStRatio.setConstant(_stancePhaseRatio);
		bias = _bias;
		startT.setConstant(_t);
		contactPast.setZero();
		phasePast.setConstant(0.5);
		statusPast = WaveStatus::STANCE_ALL;
		Tstance = period * stRatio;
		Tswing = period * (1 - stRatio);
	}

	Eigen::Vector3d GaitCtrl::calFootPos(int legID, Eigen::Vector2d vxyTargetGlobal, double dYawTarget, double phase)
	{
		this->_updateFootPoints();
		// 计算xy平面的落脚点规划
		bodyVelGlobal = bodyController->currentBalanceState.p_dot;
		bodyWGlobal = bodyController->currentBalanceState.r_dot;

		nextStep(0) = bodyVelGlobal(0) * (1 - phase) * Tswing + bodyVelGlobal(0) * Tstance / 2 + kx * (vxyTargetGlobal(0) - bodyVelGlobal(0));
		nextStep(1) = bodyVelGlobal(1) * (1 - phase) * Tswing + bodyVelGlobal(1) * Tstance / 2 + ky * (vxyTargetGlobal(1) - bodyVelGlobal(1));
		nextStep(2) = 0;

		// 计算旋转状态的落脚点叠加规划
		yaw = bodyController->currentBalanceState.r(2);
		dYaw = bodyController->currentBalanceState.r_dot(2);
		nextYaw = dYaw * (1 - phase) * Tswing + dYaw * Tstance / 2 + kyaw * (dYawTarget - dYaw);

		nextStep(0) += feetRadius(legID) * cos(yaw + feetInitAngle(legID) + nextYaw);
		nextStep(1) += feetRadius(legID) * sin(yaw + feetInitAngle(legID) + nextYaw);

		footPos = bodyController->currentBalanceState.p + nextStep;
		footPos(2) = 0.0;

		return footPos;
	}

	/*Eigen::Vector3d GaitCtrl::calFootPosW(int legID, Eigen::Vector2d vxyTargetGlobal, double dYawTarget, double phase)
	{
		Eigen::Vector3d vg = bodyController->bodyObject->est->getEstBodyVelS() - bodyController->bodyObject->est->getEstFootVelS();
		Eigen::Vector3d vsh = bodyController->bodyObject->v3_to_m3(bodyController->currentBalanceState.r_dot) * (bodyController->bodyObject->Rsb_c * initLeg[legID]);
		nextStep = 0.5 * Tstance * (vg + vsh) + kx * (vg - Eigen::Vector3d(vxyTargetGlobal(0), vxyTargetGlobal(1), 0)) + kyaw * vsh;
		nextStep += bodyController->bodyObject->v3_to_m3(0.5 * gaitHeight / 9.81 * bodyController->bodyObject->est->getEstBodyVelS()) * Eigen::Vector3d(0, 0, dYawTarget);
		yaw = bodyController->currentBalanceState.r(2);
		nextStep(0) += feetRadius(legID) * cos(yaw + feetInitAngle(legID));
		nextStep(1) += feetRadius(legID) * sin(yaw + feetInitAngle(legID));
		footPos = bodyController->currentBalanceState.p + nextStep;
		footPos(2) = 0.0;
		return footPos;
	}*/

	Eigen::Vector3d GaitCtrl::calFootPosW(int legID, Eigen::Vector2d vxyTargetGlobal, double dYawTarget, double phase, double maxStepL,const Eigen::Matrix3d& _slope)
	{
		// 计算xy平面的落脚点规划
		bodyVelGlobal = bodyController->currentBalanceState.p_dot - bodyController->bodyObject->est->getEstFootVelS();
		/*if (bodyVelGlobal.block(0, 0, 2, 1).norm() > maxStepL)
		{
			bodyVelGlobal.block(0, 0, 2, 1) = bodyVelGlobal.block(0, 0, 2, 1).normalized() * maxStepL;
		}*/
		bodyWGlobal = bodyController->currentBalanceState.r_dot;

		nextStep(0) = bodyVelGlobal(0) * (1 - phase) * Tswing + bodyVelGlobal(0) * Tstance / 2 + kx * (vxyTargetGlobal(0) - bodyController->currentBalanceState.p_dot(0));
		nextStep(1) = bodyVelGlobal(1) * (1 - phase) * Tswing + bodyVelGlobal(1) * Tstance / 2 + ky * (vxyTargetGlobal(1) - bodyController->currentBalanceState.p_dot(1));
		//nextStep(0) = bodyVelGlobal(0) * (1 - phase) * Tswing + bodyVelGlobal(0) * Tstance / 2 + kx * (vxyTargetGlobal(0) - bodyVelGlobal(0));
		//nextStep(1) = bodyVelGlobal(1) * (1 - phase) * Tswing + bodyVelGlobal(1) * Tstance / 2 + ky * (vxyTargetGlobal(1) - bodyVelGlobal(1));
		nextStep(2) = 0;

		/*if (nextStep.norm() > maxStepL)
		{
			nextStep = nextStep.normalized() * maxStepL;
		}*/

		// 计算旋转状态的落脚点叠加规划
		yaw = bodyController->currentBalanceState.r(2);
		dYaw = bodyController->currentBalanceState.r_dot(2);
		nextYaw = dYaw * (1 - phase) * Tswing + dYaw * Tstance / 2 + kyaw * (dYawTarget - dYaw);

		nextStep(0) += feetRadius(legID) * cos(yaw + feetInitAngle(legID) + nextYaw);
		nextStep(1) += feetRadius(legID) * sin(yaw + feetInitAngle(legID) + nextYaw);

		footPos = bodyController->currentBalanceState.p + nextStep;
		footPos(2) = 0.0;

		/*if (this->statusPast == WaveStatus::ADAPT)
		{
			footPos = bodyController->currentBalanceState.p + bodyController->bodyObject->Rsbh_c * initLeg[legID];
			footPos(2) = 0;
		}*/

		/*footPos = bodyController->currentBalanceState.p + _slope.transpose() * nextStep;
		footPos(2) = 0.0;*/

		return footPos;
	}

	void GaitCtrl::calcWave(Eigen::Vector4d& phase, Eigen::Vector4i& contact, WaveStatus status, double _t, Eigen::Vector4d& _estPhase, Eigen::Vector4i& _estContact)
	{
		static Vector4i lastEstContact = Vector4i::Ones();
		if (status == WaveStatus::WAVE_ALL)
		{
			passT = Eigen::Vector4d::Constant(_t) - startT;
			for (int i(0); i < 4; ++i)
			{
				// 得到总的步态周期相位
				normalT(i) = fmod(passT(i) + period - period * bias(i), period) / period;// 取余操作并归一化
				// 根据归一化的T判断应该是接触地面还是摆动
				if (normalT(i) < (1 - stRatio))
				{
					// 计算非触地过程中的相位变化0-1
					contact(i) = 0;
					phase(i) = normalT(i) / (1 - stRatio);
				}
				else
				{
					// 计算触地过程中的相位变化0-1
					contact(i) = 1;
					phase(i) = (normalT(i) + stRatio - 1) / stRatio;
				}
				// 根据被动接触检测建立相位周期
				if (_estContact(i) == 1 && lastEstContact(i) == 0)
				{
					estStRatio(i) = 1 - normalT(i);
				}
				if (normalT(i) > (1 - estStRatio(i)))
				{
					_estPhase(i) = (normalT(i) + estStRatio(i) - 1) / estStRatio(i);
				}
				else
				{
					//_estPhase(i) = normalT(i) / (1 - estStRatio(i));
					_estPhase(i) = 0;
				}
			}
		}
		else if (status == WaveStatus::SWING_ALL)
		{
			contact.setZero();
			phase.setConstant(0.5);
			_estContact.setZero();
			_estPhase.setConstant(0.5);
		}
		else if (status == WaveStatus::STANCE_ALL)
		{
			contact.setOnes();
			phase.setConstant(0.5);
			_estContact.setOnes();
			_estPhase.setConstant(0.5);
		}
		else if (status == WaveStatus::ADAPT)
		{
			passT = Eigen::Vector4d::Constant(_t) - startT;
			for (int i = 0; i < 4; i++)
			{
				if (this->recoverflag(i) == 1)
				{
					// 得到总的步态周期相位(针对每个足端的独立调整不需要考虑步态偏置)
					normalT(i) = fmod(passT(i) + 0.95 * period, period) / period;// 取余操作并归一化
					// 根据归一化的T判断应该是接触地面还是摆动
					if (normalT(i) < (1 - stRatio))
					{
						// 计算非触地过程中的相位变化0-1
						contact(i) = 0;
						phase(i) = normalT(i) / (1 - stRatio);
					}
					else
					{
						// 计算触地过程中的相位变化0-1
						contact(i) = 1;
						phase(i) = (normalT(i) + stRatio - 1) / stRatio;
						this->recoverflag(i) = 0;
					}
				}
				else
				{
					contact(i) = 1;
					phase(i) = 0.5;
				}
				
				//if (passT(i) > 0.95 * period)
				//{
				//	// 完整周期结束之后，自动关闭恢复相位
				//	this->recoverflag(i) = 0;
				//}
			}
		}
		lastEstContact = _estContact;
	}

	void GaitCtrl::calcContactPhase(WaveStatus status, double _t, Eigen::Vector4d& _estPhase, Eigen::Vector4i& _estContact)
	{

		calcWave(phase, contact, status, _t, _estPhase, _estContact);

		if (status != statusPast)
		{
			if (switchStatus.sum() == 0)
			{
				switchStatus.setOnes();
			}
			calcWave(phasePast, contactPast, statusPast, _t, _estPhase, _estContact);
			// 两种情况，分别是从完全站立到全摆动，以及全摆动到完全站立
			if ((status == WaveStatus::STANCE_ALL) && (statusPast == WaveStatus::SWING_ALL))
			{
				contactPast.setOnes();
			}
			else if ((status == WaveStatus::SWING_ALL) && (statusPast == WaveStatus::STANCE_ALL))
			{
				contactPast.setZero();
			}
			else if ((status == WaveStatus::ADAPT) && (statusPast == WaveStatus::SWING_ALL))
			{
				contactPast.setOnes();
			}
			else if ((status == WaveStatus::SWING_ALL) && (statusPast == WaveStatus::ADAPT))
			{
				contactPast.setOnes();
			}
			else if ((status == WaveStatus::ADAPT) && (statusPast == WaveStatus::STANCE_ALL))
			{
				contactPast.setOnes();
			}
			else if ((status == WaveStatus::SWING_ALL) && (statusPast == WaveStatus::ADAPT))
			{
				contactPast.setOnes();
			}
		}

		// 如果切换状态为允许
		if (switchStatus.sum() != 0)
		{
			for (int i = 0; i < 4; ++i)
			{
				if (contact(i) == contactPast(i))
				{
					// 切换完成
					switchStatus(i) = 0;
				}
				else
				{
					contact(i) = contactPast(i);
					phase(i) = phasePast(i);
				}
			}
			if (switchStatus.sum() == 0)
			{
				statusPast = status;
			}
		}

		*gaitPhase = phase;
		*gaitContact = contact;
	}

	double GaitCtrl::cycloidXYPosition(double start, double end, double phase) {
		double phasePI = 2 * M_PI * phase;
		return (end - start) * (phasePI - sin(phasePI)) / (2 * M_PI) + start;
	}

	double GaitCtrl::cycloidXYVelocity(double start, double end, double phase) {
		double phasePI = 2 * M_PI * phase;
		return (end - start) * (1 - cos(phasePI)) / Tswing;
	}

	double GaitCtrl::cycloidZPosition(double start, double h, double phase) {
		double phasePI = 2 * M_PI * phase;
		return h * (1 - cos(phasePI)) / 2 + start;
	}

	double GaitCtrl::cycloidZVelocity(double h, double phase) {
		double phasePI = 2 * M_PI * phase;
		return h * M_PI * sin(phasePI) / Tswing;
	}

	void GaitCtrl::setGait(Eigen::Vector2d _VxyTargetGlobal, double _dYawTarget, double _gaitHeight)
	{
		this->VxyTarget = _VxyTargetGlobal;
		this->dYawTarget = _dYawTarget;
		this->gaitHeight = _gaitHeight;
	}

	void GaitCtrl::run(Eigen::Matrix<double, 3, 4>& _feetPos, Eigen::Matrix<double, 3, 4>& _feetVel, double _maxStepL, const Eigen::Matrix3d& _slope)
	{
		if (is_firstRun) {
			//this->startT.setConstant(_t);
			//startP = bodyController->bodyObject->est->getEstFeetPosS();
			startP = bodyController->bodyObject->getFKFeetPos();
			is_firstRun = false;
		}

		for (int i(0); i < 4; ++i) {
			if ((*gaitContact)(i) == 1) {
				//if ((*gaitPhase)(i) < 0.5) {
					//startP.col(i) = bodyController->bodyObject->est->getEstFeetPosS(i);
					startP.col(i) = bodyController->bodyObject->getFKFeetPos(i);
				//}
				_feetPos.col(i) = startP.col(i);
				_feetVel.col(i).setZero();
			}
			else {
				// 自适应的目标是让足端回到当前中性立足点，要消除速度的影响
				/*if (this->statusPast == WaveStatus::ADAPT)
				{
					VxyTarget = bodyController->currentBalanceState.p_dot.segment(0, 2);
					dYawTarget = bodyController->currentBalanceState.r_dot(2);
				}*/
				endP.col(i) = calFootPosW(i, VxyTarget, dYawTarget, (*gaitPhase)(i), _maxStepL, _slope);

				_feetPos.col(i) = getFootPos(i);
				_feetVel.col(i) = getFootVel(i);
			}
		}
		pastP = _feetPos;
		gaitPhasePast = *gaitPhase;
	}

	Eigen::Vector3d GaitCtrl::getFootPos(int i)
	{
		Eigen::Vector3d footPos;

		footPos(0) = cycloidXYPosition(startP.col(i)(0), endP.col(i)(0), (*gaitPhase)(i));
		footPos(1) = cycloidXYPosition(startP.col(i)(1), endP.col(i)(1), (*gaitPhase)(i));
		footPos(2) = cycloidZPosition(startP.col(i)(2), gaitHeight, (*gaitPhase)(i));

		return footPos;
	}

	Eigen::Vector3d GaitCtrl::getFootVel(int i)
	{
		Eigen::Vector3d footVel;

		footVel(0) = cycloidXYVelocity(startP.col(i)(0), endP.col(i)(0), (*gaitPhase)(i));
		footVel(1) = cycloidXYVelocity(startP.col(i)(1), endP.col(i)(1), (*gaitPhase)(i));
		footVel(2) = cycloidZVelocity(gaitHeight, (*gaitPhase)(i));

		return footVel;
	}

	void GaitCtrl::restart()
	{
		this->is_firstRun = true;
		this->VxyTarget.setZero();
	}

	void GaitCtrl::footPointMarking()
	{
		Eigen::Matrix<double, 3, 4> cur = this->bodyController->bodyObject->getFKFeetPosB();
		Eigen::Matrix<double, 3, 4> tar = this->bodyController->bodyObject->initLegsXYPosition;
		Eigen::Matrix<double, 3, 4> err = tar - cur;
		double px = 0.2,py=0.1;
		for (int i = 0; i < 4; i++)
		{
			// eth marking
			this->fpmark(i) = 1 - sqrt(pow(err(0, i) / px, 2) + pow(err(1, i) / py, 2));
			// other marking
			//this->fpmark(i) = exp(-err.block(0, i, 1, 1).squaredNorm() / (2 * px * px) - err.block(1, i, 1, 1).squaredNorm() / (2 * py * py));
		}
		std::cout << "mark: " << this->fpmark.transpose() << std::endl;
	}

	void GaitCtrl::footStateTransform(double _t)
	{
		double threshold = 0.;
		bool required_transform = false;
		double low = this->fpmark(0);
		uint8_t low_idx = 0;
		std::vector<int> lowMark;
		lowMark.clear();
		// 找出打分最低的那个
		for (int i = 0; i < 4; i++)
		{
			if (!this->recoverflag(i))
			{
				if (this->fpmark(i) < low)
				{
					low = this->fpmark(i);
					low_idx = i;
				}
				if (this->fpmark(i) < threshold)
				{
					if (low_idx == i)
					{
						// 只需要优先处理打分最低的即可
						lowMark.insert(lowMark.begin(), i);
					}
					else
					{
						lowMark.push_back(i);
					}
				}
			}
		}
		// 若打分最低的超过了阈值，则需要调整
		if (this->statusPast == WaveStatus::ADAPT)
		{
			if (low < threshold)
			{
				for (int i = 0; i < lowMark.size(); i++)
				{
					// 需要保证相邻的两足端不会同时调整
					switch (lowMark.at(i)) {
					case 0:
						if (!this->recoverflag(1) && !this->recoverflag(2))
						{
							this->recoverflag(0) = 1;
						}
						break;
					case 1:
						if (!this->recoverflag(0) && !this->recoverflag(3))
						{
							this->recoverflag(1) = 1;
						}
						break;
					case 2:
						if (!this->recoverflag(0) && !this->recoverflag(3))
						{
							this->recoverflag(2) = 1;
						}
						break;
						break;
					case 3:
						if (!this->recoverflag(1) && !this->recoverflag(2))
						{
							this->recoverflag(3) = 1;
						}
						break;
					default:
						break;
					}
				}
			}
		}
		else
		{
			this->recoverflag.setZero();
		}
		
		// 如果开启了调整需要同步更新起始时间
		for (int i = 0; i < 4; i++)
		{
			if (this->recoverflag(i) != this->recoverflaglast(i))
			{
				this->startT(i) = _t;
			}
		}

		this->recoverflaglast = this->recoverflag;
		std::cout << "flag: " << this->recoverflag.transpose() << std::endl;
	}
}
