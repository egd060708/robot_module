#pragma once
#include "kelmanFilter.h"

class EstBase {
protected:
	Eigen::MatrixXd _H;
	Eigen::MatrixXd _A;
	Eigen::MatrixXd _B;
	Eigen::MatrixXd _R;
	Eigen::MatrixXd _RInit;
	Eigen::VectorXd _Rdig;
	Eigen::MatrixXd _Q;
	Eigen::MatrixXd _QInit;
	Eigen::VectorXd _Qdig;
	Eigen::MatrixXd _Cu;
	Eigen::MatrixXd _P;

	Eigen::Matrix4d _Tsb;
public:
	Eigen::VectorXd estimatorOut;
	Eigen::VectorXd estimatorState;
	inline virtual void estimatorRun(const Eigen::MatrixXd& _u,const Eigen::MatrixXd& _y,const Eigen::Vector4i& _contact,const Eigen::Vector4d& _phase) = 0;
	inline virtual Eigen::Vector3d getEstFeetPosS(int id) = 0;
	inline virtual Eigen::Vector3d getEstFeetVelS(int id) = 0;
	inline virtual Eigen::Matrix<double, 3, 4> getEstFeetPosS() = 0;
	inline virtual Eigen::Matrix<double, 3, 4> getEstFeetVelS() = 0;
	inline virtual Eigen::Vector3d getEstBodyPosS() = 0;
	inline virtual Eigen::Vector3d getEstBodyVelS() = 0;
	inline void updateTsb(const Eigen::MatrixXd& _T)
	{
		_Tsb = _T;
	}
};

class QpEst : public EstBase {
private:
	kelmanFilter<18, 3, 28> estimator;
	double _largeVariance = 100;// 大的协方差
	double _trust;// 对于腿部是否触地的置信度
	template<typename T>
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
	}
public:
	QpEst(double _dt) {
		_H.resize(28, 18);
		_A.resize(18, 18);
		_B.resize(18, 3);
		_R.resize(28, 28);
		_RInit.resize(28, 28);
		_Rdig.resize(28);
		_Q.resize(18, 18);
		_QInit.resize(18, 18);
		_Qdig.resize(18);
		_Cu.resize(3, 3);
		_P.resize(18, 18);
		estimatorOut.resize(28);
		estimatorState.resize(18);

		_H.setZero();
		_A.setZero();
		_B.setZero();
		_R.setZero();
		_RInit.setZero();
		_Q.setZero();
		_QInit.setZero();
		_Qdig.setZero();
		_Cu.setZero();
		_P.setIdentity();
		estimatorOut.setZero();
		estimatorState.setZero();

		// 初始化状态估计器
		// 状态转移矩阵
		_A.block(0, 3, 3, 3) = Eigen::Matrix<double, 3, 3>::Identity();
		// 输入矩阵
		_B.block(3, 0, 3, 3) = Eigen::Matrix<double, 3, 3>::Identity();
		// 输出矩阵
		_H.block(0, 0, 3, 3) = -Eigen::Matrix<double, 3, 3>::Identity();
		_H.block(3, 0, 3, 3) = -Eigen::Matrix<double, 3, 3>::Identity();
		_H.block(6, 0, 3, 3) = -Eigen::Matrix<double, 3, 3>::Identity();
		_H.block(9, 0, 3, 3) = -Eigen::Matrix<double, 3, 3>::Identity();
		_H.block(12, 3, 3, 3) = -Eigen::Matrix<double, 3, 3>::Identity();
		_H.block(15, 3, 3, 3) = -Eigen::Matrix<double, 3, 3>::Identity();
		_H.block(18, 3, 3, 3) = -Eigen::Matrix<double, 3, 3>::Identity();
		_H.block(21, 3, 3, 3) = -Eigen::Matrix<double, 3, 3>::Identity();
		_H.block(0, 6, 12, 12) = Eigen::Matrix<double, 12, 12>::Identity();
		_H(24, 8) = 1;
		_H(25, 11) = 1;
		_H(26, 14) = 1;
		_H(27, 17) = 1;
		// 初始化状态空间方程
		estimator.setFunc(_A, _B, _H, _dt);
		// 测量噪声协方差
		_RInit << 0.008, 0.012, -0.000, -0.009, 0.012, 0.000, 0.009, -0.009, -0.000, -0.009, -0.009, 0.000, -0.000, 0.000, -0.000, 0.000, -0.000, -0.001, -0.002, 0.000, -0.000, -0.003, -0.000, -0.001, 0.000, 0.000, 0.000, 0.000,
			0.012, 0.019, -0.001, -0.014, 0.018, -0.000, 0.014, -0.013, -0.000, -0.014, -0.014, 0.001, -0.001, 0.001, -0.001, 0.000, 0.000, -0.001, -0.003, 0.000, -0.001, -0.004, -0.000, -0.001, 0.000, 0.000, 0.000, 0.000,
			-0.000, -0.001, 0.001, 0.001, -0.001, 0.000, -0.000, 0.000, -0.000, 0.001, 0.000, -0.000, 0.000, -0.000, 0.000, 0.000, -0.000, -0.000, 0.000, -0.000, -0.000, -0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
			-0.009, -0.014, 0.001, 0.010, -0.013, 0.000, -0.010, 0.010, 0.000, 0.010, 0.010, -0.000, 0.001, 0.000, 0.000, 0.001, -0.000, 0.001, 0.002, -0.000, 0.000, 0.003, 0.000, 0.001, 0.000, 0.000, 0.000, 0.000,
			0.012, 0.018, -0.001, -0.013, 0.018, -0.000, 0.013, -0.013, -0.000, -0.013, -0.013, 0.001, -0.001, 0.000, -0.001, 0.000, 0.001, -0.001, -0.003, 0.000, -0.001, -0.004, -0.000, -0.001, 0.000, 0.000, 0.000, 0.000,
			0.000, -0.000, 0.000, 0.000, -0.000, 0.001, 0.000, 0.000, -0.000, 0.000, 0.000, -0.000, -0.000, 0.000, -0.000, 0.000, 0.000, 0.000, -0.000, -0.000, -0.000, -0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
			0.009, 0.014, -0.000, -0.010, 0.013, 0.000, 0.010, -0.010, -0.000, -0.010, -0.010, 0.000, -0.001, 0.000, -0.001, 0.000, -0.000, -0.001, -0.001, 0.000, -0.000, -0.003, -0.000, -0.001, 0.000, 0.000, 0.000, 0.000,
			-0.009, -0.013, 0.000, 0.010, -0.013, 0.000, -0.010, 0.009, 0.000, 0.010, 0.010, -0.000, 0.001, -0.000, 0.000, -0.000, 0.000, 0.001, 0.002, 0.000, 0.000, 0.003, 0.000, 0.001, 0.000, 0.000, 0.000, 0.000,
			-0.000, -0.000, -0.000, 0.000, -0.000, -0.000, -0.000, 0.000, 0.001, 0.000, 0.000, 0.000, 0.000, -0.000, 0.000, -0.000, 0.000, -0.000, 0.000, -0.000, 0.000, 0.000, -0.000, -0.000, 0.000, 0.000, 0.000, 0.000,
			-0.009, -0.014, 0.001, 0.010, -0.013, 0.000, -0.010, 0.010, 0.000, 0.010, 0.010, -0.000, 0.001, 0.000, 0.000, -0.000, -0.000, 0.001, 0.002, -0.000, 0.000, 0.003, 0.000, 0.001, 0.000, 0.000, 0.000, 0.000,
			-0.009, -0.014, 0.000, 0.010, -0.013, 0.000, -0.010, 0.010, 0.000, 0.010, 0.010, -0.000, 0.001, -0.000, 0.000, -0.000, 0.000, 0.001, 0.002, -0.000, 0.000, 0.003, 0.001, 0.001, 0.000, 0.000, 0.000, 0.000,
			0.000, 0.001, -0.000, -0.000, 0.001, -0.000, 0.000, -0.000, 0.000, -0.000, -0.000, 0.001, 0.000, -0.000, -0.000, -0.000, 0.000, 0.000, -0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
			-0.000, -0.001, 0.000, 0.001, -0.001, -0.000, -0.001, 0.001, 0.000, 0.001, 0.001, 0.000, 1.708, 0.048, 0.784, 0.062, 0.042, 0.053, 0.077, 0.001, -0.061, 0.046, -0.019, -0.029, 0.000, 0.000, 0.000, 0.000,
			0.000, 0.001, -0.000, 0.000, 0.000, 0.000, 0.000, -0.000, -0.000, 0.000, -0.000, -0.000, 0.048, 5.001, -1.631, -0.036, 0.144, 0.040, 0.036, 0.016, -0.051, -0.067, -0.024, -0.005, 0.000, 0.000, 0.000, 0.000,
			-0.000, -0.001, 0.000, 0.000, -0.001, -0.000, -0.001, 0.000, 0.000, 0.000, 0.000, -0.000, 0.784, -1.631, 1.242, 0.057, -0.037, 0.018, 0.034, -0.017, -0.015, 0.058, -0.021, -0.029, 0.000, 0.000, 0.000, 0.000,
			0.000, 0.000, 0.000, 0.001, 0.000, 0.000, 0.000, -0.000, -0.000, -0.000, -0.000, -0.000, 0.062, -0.036, 0.057, 6.228, -0.014, 0.932, 0.059, 0.053, -0.069, 0.148, 0.015, -0.031, 0.000, 0.000, 0.000, 0.000,
			-0.000, 0.000, -0.000, -0.000, 0.001, 0.000, -0.000, 0.000, 0.000, -0.000, 0.000, 0.000, 0.042, 0.144, -0.037, -0.014, 3.011, 0.986, 0.076, 0.030, -0.052, -0.027, 0.057, 0.051, 0.000, 0.000, 0.000, 0.000,
			-0.001, -0.001, -0.000, 0.001, -0.001, 0.000, -0.001, 0.001, -0.000, 0.001, 0.001, 0.000, 0.053, 0.040, 0.018, 0.932, 0.986, 0.885, 0.090, 0.044, -0.055, 0.057, 0.051, -0.003, 0.000, 0.000, 0.000, 0.000,
			-0.002, -0.003, 0.000, 0.002, -0.003, -0.000, -0.001, 0.002, 0.000, 0.002, 0.002, -0.000, 0.077, 0.036, 0.034, 0.059, 0.076, 0.090, 6.230, 0.139, 0.763, 0.013, -0.019, -0.024, 0.000, 0.000, 0.000, 0.000,
			0.000, 0.000, -0.000, -0.000, 0.000, -0.000, 0.000, 0.000, -0.000, -0.000, -0.000, 0.000, 0.001, 0.016, -0.017, 0.053, 0.030, 0.044, 0.139, 3.130, -1.128, -0.010, 0.131, 0.018, 0.000, 0.000, 0.000, 0.000,
			-0.000, -0.001, -0.000, 0.000, -0.001, -0.000, -0.000, 0.000, 0.000, 0.000, 0.000, 0.000, -0.061, -0.051, -0.015, -0.069, -0.052, -0.055, 0.763, -1.128, 0.866, -0.022, -0.053, 0.007, 0.000, 0.000, 0.000, 0.000,
			-0.003, -0.004, -0.000, 0.003, -0.004, -0.000, -0.003, 0.003, 0.000, 0.003, 0.003, 0.000, 0.046, -0.067, 0.058, 0.148, -0.027, 0.057, 0.013, -0.010, -0.022, 2.437, -0.102, 0.938, 0.000, 0.000, 0.000, 0.000,
			-0.000, -0.000, 0.000, 0.000, -0.000, 0.000, -0.000, 0.000, -0.000, 0.000, 0.001, 0.000, -0.019, -0.024, -0.021, 0.015, 0.057, 0.051, -0.019, 0.131, -0.053, -0.102, 4.944, 1.724, 0.000, 0.000, 0.000, 0.000,
			-0.001, -0.001, 0.000, 0.001, -0.001, 0.000, -0.001, 0.001, -0.000, 0.001, 0.001, 0.000, -0.029, -0.005, -0.029, -0.031, 0.051, -0.003, -0.024, 0.018, 0.007, 0.938, 1.724, 1.569, 0.000, 0.000, 0.000, 0.000,
			0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 1.0, 0.000, 0.000, 0.000,
			0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 1.0, 0.000, 0.000,
			0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 1.0, 0.000,
			0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 1.0;
		// 过程噪声协方差
		for (int i(0); i < _Qdig.rows(); ++i)
		{
			if (i < 3)
			{
				_Qdig(i) = 0.0003;
			}
			else if (i < 6)
			{
				_Qdig(i) = 0.0003;
			}
			else
			{
				_Qdig(i) = 0.01;
			}
		}

		_Cu << 268.573, -43.819, -147.211,
			-43.819, 92.949, 58.082,
			-147.211, 58.082, 302.120;
		// 过程协方差
		_QInit = _Qdig.asDiagonal();
		_QInit += _B * _Cu * _B.transpose();// 把加速度计的协方差矩阵嵌入过程协方差矩阵中

		estimator.setConv(_QInit, _RInit, _largeVariance * _P);// 初始化一个比较大的预测协方差矩阵
	}

	inline void estimatorRun(const Eigen::MatrixXd& _u,const Eigen::MatrixXd& _y, const Eigen::Vector4i& _contact, const Eigen::Vector4d& _phase) override
	{
		_R = _RInit;
		_Q = _QInit;
		for (int i = 0; i < 4; i++)
		{
			if (_contact(i) == 0)
			{
				_Q.block(6 + 3 * i, 6 + 3 * i, 3, 3) = _largeVariance * Eigen::Matrix<double, 3, 3>::Identity();
				_R.block(12 + 3 * i, 12 + 3 * i, 3, 3) = _largeVariance * Eigen::Matrix<double, 3, 3>::Identity();
				_R(24 + i, 24 + i) = _largeVariance;
			}
			else
			{
				_trust = windowFunc(_phase(i), 0.2);
				_Q.block(6 + 3 * i, 6 + 3 * i, 3, 3) = (1 + (1 - _trust) * _largeVariance) * _QInit.block(6 + 3 * i, 6 + 3 * i, 3, 3);
				//std::cout << "trustM: " << _Q.block(6 + 3 * i, 6 + 3 * i, 3, 3) << std::endl;
				_R.block(12 + 3 * i, 12 + 3 * i, 3, 3) = (1 + (1 - _trust) * _largeVariance) * _RInit.block(12 + 3 * i, 12 + 3 * i, 3, 3);
				_R(24 + i, 24 + i) = (1 + (1 - _trust) * _largeVariance) * _RInit(24 + i, 24 + i);
			}
		}
		// 更新协方差矩阵
		estimator.updateConv(_Q, _R);
		// 卡尔曼滤波执行
		estimator.f(_u, _y);
		// 估计输出
		estimatorOut = estimator.getOut();
		// 估计状态
		estimatorState = estimator.getState();
	}

	inline Eigen::Vector3d getEstFeetPosS(int id) override
	{
		Eigen::Vector3d out;
		out = estimatorState.block<3, 1>(6 + 3 * id, 0);
		return out;
	}

	inline Eigen::Vector3d getEstFeetVelS(int id) override
	{
		Eigen::Vector3d out;
		out = estimatorOut.block<3, 1>(12 + 3 * id, 0);
		return out;
	}

	inline Eigen::Matrix<double, 3, 4> getEstFeetPosS() override
	{
		Eigen::Matrix<double, 3, 4> out;
		for (int i = 0; i < 4; i++)
		{
			out.block<3, 1>(0, i) = estimatorState.block<3, 1>(6 + 3 * i, 0);
		}
		return out;
	}

	inline Eigen::Matrix<double, 3, 4> getEstFeetVelS() override
	{
		Eigen::Matrix<double, 3, 4> out;
		for (int i = 0; i < 4; i++)
		{
			out.block<3, 1>(0, i) = estimatorOut.block<3, 1>(12 + 3 * i, 0);
		}
		return out;
	}

	inline Eigen::Vector3d getEstBodyPosS() override
	{
		return estimatorState.block<3, 1>(0, 0);
	}
	
	inline Eigen::Vector3d getEstBodyVelS() override
	{
		return estimatorState.block<3, 1>(3, 0);
	}
};

class QpwEst : public EstBase {
private:
	kelmanFilter<24, 3, 44> estimator;
	double _largeVariance = 100;// 大的协方差
	double _ctTrust = 0;// 对于腿部是否触地的置信度
	double _czTrust = 0;// 对于腿部参考高度的置信度
	double wp=0.0003, wpd=0.0003, wpw=0.0003, wpwd=0.0003, wpcp=0.01;// 过程协方差参数
	double vpcp=0.01, vpcpkd=0.01, vpcpwd=0.01, wz = 1.0, cpz = 1.0;// 测量协方差参数
	/* 非线性函数 */
	double sigmoid(double _x) {
		return 1 / (1 + exp(-_x));
	}
	/* 触地状态信任窗口函数 */
	double ctWin(double _phase, double _w=0.1)
	{
		return 0.5 * (sigmoid(4 * _phase / _w - 2) + sigmoid(4 * (1 - _phase) / _w - 2));
	}
	/* 高度设定信任窗口函数 */
	double czWin(double _z, double _k1 = 50, double _k2 = 10)
	{
		if (_z > 0)
		{
			return exp(-_k1 * _z * _z);
		}
		else
		{
			return exp(-_k2 * _z * _z);
		}
	}
	template<typename T>
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
	}
public:
	QpwEst(double _dt) {
		_H.resize(44, 24);
		_A.resize(24, 24);
		_B.resize(24, 3);
		_R.resize(44, 44);
		_RInit.resize(44, 44);
		_Rdig.resize(44);
		_Q.resize(24, 24);
		_QInit.resize(24, 24);
		_Qdig.resize(24);
		_Cu.resize(3, 3);
		_P.resize(24, 24);
		estimatorOut.resize(44);
		estimatorState.resize(24);

		_H.setZero();
		_A.setZero();
		_B.setZero();
		_R.setZero();
		_RInit.setZero();
		_Q.setZero();
		_QInit.setZero();
		_Qdig.setZero();
		_Cu.setZero();
		_P.setIdentity();
		estimatorOut.setZero();
		estimatorState.setZero();

		// 初始化观测器参数
		_A.block(0, 3, 3, 3) = Eigen::Matrix3d::Identity();
		_A.block(6, 9, 3, 3) = Eigen::Matrix3d::Identity();
		_B.block(3, 0, 3, 3) = Eigen::Matrix3d::Identity();

		for (int i = 0; i < 4; i++)
		{
			_H.block(3 * i, 0, 3, 3) = -Eigen::Matrix3d::Identity();
			_H.block(3 * i, 6, 3, 3) = Eigen::Matrix3d::Identity();
			_H.block(12 + 3 * i, 3, 3, 3) = -Eigen::Matrix3d::Identity();
			_H.block(12 + 3 * i, 9, 3, 3) = Eigen::Matrix3d::Identity();
			_H.block(24 + 3 * i, 9, 3, 3) = Eigen::Matrix3d::Identity();
			_H.block(36 + i, 6, 1, 3) = Eigen::Vector3d(0, 0, 1).transpose();
			_H.block(40 + i, 12 + 3 * i, 1, 3) = Eigen::Vector3d(0, 0, 1).transpose();
		}
		_H.block(0, 12, 12, 12) = Eigen::Matrix<double, 12, 12>::Identity();

		estimator.setFunc(_A, _B, _H, _dt);

		_Qdig.block(0, 0, 3, 1).setConstant(wp);
		_Qdig.block(3, 0, 3, 1).setConstant(wpd);
		_Qdig.block(6, 0, 3, 1).setConstant(wpw);
		_Qdig.block(9, 0, 3, 1).setConstant(wpwd);
		_Qdig.block(12, 0, 12, 1).setConstant(wpcp);

		_Rdig.block(0, 0, 12, 1).setConstant(vpcp);
		_Rdig.block(12, 0, 12, 1).setConstant(vpcpkd);
		_Rdig.block(24, 0, 12, 1).setConstant(vpcpwd);
		_Rdig.block(36, 0, 4, 1).setConstant(wz);
		_Rdig.block(40, 0, 4, 1).setConstant(cpz);

		_Cu << 268.573, -43.819, -147.211,
			-43.819, 92.949, 58.082,
			-147.211, 58.082, 302.120;
		// 过程协方差
		_QInit = _Qdig.asDiagonal();
		_QInit += _B * _Cu * _B.transpose();// 把加速度计的协方差矩阵嵌入过程协方差矩阵中
		_RInit = _Rdig.asDiagonal();

		estimator.setConv(_QInit, _RInit, _largeVariance * _P);// 初始化一个比较大的预测协方差矩阵
	}

	inline void estimatorRun(const Eigen::MatrixXd& _u, const Eigen::MatrixXd& _y, const Eigen::Vector4i& _contact, const Eigen::Vector4d& _phase) override
	{
		_R = _RInit;
		_Q = _QInit;
		for (int i = 0; i < 4; i++)
		{
			if (_contact(i) == 0)
			{
				_Q.block(12 + 3 * i, 12 + 3 * i, 3, 3) = _largeVariance * _QInit.block(12 + 3 * i, 12 + 3 * i, 3, 3);
				_R.block(3 * i, 3 * i, 3, 3) = _largeVariance * _RInit.block(3 * i, 3 * i, 3, 3);
				_R.block(12 + 3 * i, 12 + 3 * i, 3, 3) = _largeVariance * _RInit.block(12 + 3 * i, 12 + 3 * i, 3, 3);
				_R.block(24 + 3 * i, 24 + 3 * i, 3, 3) = _largeVariance * _RInit.block(24 + 3 * i, 24 + 3 * i, 3, 3);
				_R(36 + i, 36 + i) = _largeVariance * _RInit(36 + i, 36 + i);
				_R(40 + i, 40 + i) = _largeVariance * _RInit(40 + i, 40 + i);
			}
			else
			{
				_ctTrust = ctWin(_phase(i), 0.1);// 针对触地相位做信任度估计
				/*_ctTrust = windowFunc(_phase(i), 0.2);*/
				//_czTrust = czWin(_y(3 * i + 2));// 针对当前设置高度做信任度估计
				_czTrust = czWin(0);// 针对当前设置高度做信任度估计(暂且设置为0)
				Eigen::Vector3d _trustM(1, 1, _czTrust);// 获取组合置信度
				_trustM = _ctTrust * _trustM;
				_trustM = Eigen::Vector3d::Ones() + _largeVariance * (Eigen::Vector3d::Ones() - _trustM);
				_Q.block(12 + 3 * i, 12 + 3 * i, 3, 3) = _trustM.asDiagonal() * _QInit.block(12 + 3 * i, 12 + 3 * i, 3, 3);
				//std::cout << "trustM: " << _Q.block(12 + 3 * i, 12 + 3 * i, 3, 3) << std::endl;
				_R.block(3 * i, 3 * i, 3, 3) = _trustM.asDiagonal() * _RInit.block(3 * i, 3 * i, 3, 3);
				_R.block(12 + 3 * i, 12 + 3 * i, 3, 3) = _trustM.asDiagonal() * _RInit.block(12 + 3 * i, 12 + 3 * i, 3, 3);
				_R.block(24 + 3 * i, 24 + 3 * i, 3, 3) = _trustM.asDiagonal() * _RInit.block(24 + 3 * i, 24 + 3 * i, 3, 3);
				_R(36 + i, 36 + i) = _trustM(2) * _RInit(36 + i, 36 + i);
				_R(40 + i, 40 + i) = _trustM(2) * _RInit(40 + i, 40 + i);
			}
		}
		/*std::cout << "_R: \n" << _R << std::endl;
		std::cout << "_Q: \n" << _Q << std::endl;*/
		// 更新协方差矩阵
		estimator.updateConv(_Q, _R);
		// 卡尔曼滤波执行
		estimator.f(_u, _y);
		// 估计输出
		estimatorOut = estimator.getOut();
		// 估计状态
		estimatorState = estimator.getState();
	}

	inline Eigen::Vector3d getEstFeetPosS(int id) override
	{
		////Eigen::Vector3d out;
		//Eigen::Vector4d trans(0,0,0,1);
		//trans.block(0, 0, 3, 1) = estimatorOut.block<3, 1>(3 * id, 0);
		//trans = this->_Tsb * trans;
		////out = trans.block(0, 0, 3, 1);
		//return trans.block(0, 0, 3, 1);

		Eigen::Vector3d out;
		out = estimatorState.block<3, 1>(12 + 3 * id, 0);
		return out;
	}

	inline Eigen::Vector3d getEstFeetVelS(int id) override
	{
		Eigen::Vector3d out;
		//out = this->_Tsb.block(0,0,3,3) * (estimatorOut.block<3, 1>(12 + 3 * id, 0)+estimatorOut.block<3,1>(24 + 3 * id, 0));
		out = estimatorOut.block<3, 1>(12 + 3 * id, 0) + estimatorOut.block<3, 1>(24 + 3 * id, 0);
		//out = estimatorOut.block<3, 1>(12 + 3 * id, 0);
		return out;
	}

	inline Eigen::Matrix<double, 3, 4> getEstFeetPosS() override
	{
		Eigen::Matrix<double, 3, 4> out;
		
		for (int i = 0; i < 4; i++)
		{
			/*Eigen::Vector4d trans(0,0,0,1);
			trans.block(0, 0, 3, 1) = estimatorOut.block<3, 1>(3 * i, 0);
			trans = this->_Tsb * trans;
			out.block<3, 1>(0, i) = trans.block(0,0,3,1);*/

			out.block(0, 0, 3, 1) = estimatorState.block<3, 1>(12 + 3 * i, 0);
		}
		return out;
	}

	inline Eigen::Matrix<double, 3, 4> getEstFeetVelS() override
	{
		Eigen::Matrix<double, 3, 4> out;
		for (int i = 0; i < 4; i++)
		{
			//out.block<3, 1>(0, i) = this->_Tsb.block(0, 0, 3, 3) * (estimatorOut.block<3, 1>(12 + 3 * i, 0) + estimatorOut.block<3, 1>(24 + 3 * i, 0));
			out.block<3, 1>(0, i) = estimatorOut.block<3, 1>(12 + 3 * i, 0) + estimatorOut.block<3,1>(24 + 3 * i, 0);
			//out.block<3, 1>(0, i) = estimatorOut.block<3, 1>(12 + 3 * i, 0);
		}
		return out;
	}

	inline Eigen::Vector3d getEstBodyPosS() override
	{
		return estimatorState.block<3, 1>(0, 0);
	}

	inline Eigen::Vector3d getEstBodyVelS() override
	{
		return estimatorState.block<3, 1>(3, 0);
	}

	inline Eigen::Vector3d getEstFootPosS()
	{
		return estimatorState.block<3, 1>(6, 0);
	}

	inline Eigen::Vector3d getEstFootVelS()
	{
		return estimatorState.block<3, 1>(9, 0);
	}
};