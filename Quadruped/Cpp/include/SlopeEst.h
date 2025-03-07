#pragma once

#include <Eigen/Dense>


namespace Quadruped {
	class SlopeEst {
	public:
		Eigen::MatrixXd points;// 用于估计平面的点
		Eigen::MatrixXd estFaces;// 计算得到的原始平面
		Eigen::MatrixXd estFacesNormal;// 计算得到的原始平面法向量
		Eigen::Vector3d estNormal;// 估计平面得到的归一化法向量
		Eigen::Vector4i usedPoints;// 计算时使用的点使能序列

		/* 根据情况计算估计平面 */
		void _calFaces() {
			int points_num = 0;
			for (int i=0;i<4;i++)
			{
				if (usedPoints(i) == 1)
				{
					points_num++;
				}
			}

			if (points_num < 3)
			{
				// 如果不满足估计平面的条件，则直接定义为法向竖直朝上的平面
				estFaces.setZero();
				estFaces.row(2).setConstant(1);
			}
			else if (points_num == 3)
			{
				Eigen::Matrix<double, 4, 1> one_face;
				Eigen::Matrix<double, 3, 3> three_points;
				for (int i = 0; i < 4; i++)
				{
					if (usedPoints(i) == 1)
					{
						three_points.col(i) = this->points.col(i);
					}
				}
				one_face = _calOneFace(three_points);
				this->estFaces.setZero();
				// 保证法向量方向向上
				if (one_face(2) < 0)
				{
					this->estFaces.col(0) = -one_face;
				}
				else
				{
					this->estFaces.col(0) = one_face;
				}
			}
			else if (points_num == 4)
			{
				Eigen::Matrix<double, 4, 1> one_face;
				Eigen::Matrix<double, 3, 3> three_points;
				// 注意为了使求得的平面法向方向相同，此处的多点索引旋转方向应该一致，此处为逆时针
				Eigen::Matrix<int, 3, 4> idx;
				idx.col(0) << 0, 2, 1;
				idx.col(1) << 0, 3, 1;
				idx.col(2) << 0, 2, 3;
				idx.col(3) << 1, 2, 3;
				for (int i = 0; i < 4; i++)
				{
					three_points.col(0) = this->points.col(idx(0, i));
					three_points.col(1) = this->points.col(idx(1, i));
					three_points.col(2) = this->points.col(idx(2, i));
					one_face = _calOneFace(three_points);
					// 保证法向量方向为向上
					if (one_face(2) < 0)
					{
						this->estFaces.col(i) = -one_face;
					}
					else
					{
						this->estFaces.col(i) = one_face;
					}
				}
			}
			// 计算得到四个平面之后把法向量提出
			for (int i = 0; i < 4; i++)
			{
				this->estFacesNormal = this->estFaces.block(0, 0, 3, 4);
			}
		}

		/* 通过三个点计算三个平面的参数 */
		Eigen::Matrix<double, 4, 1> _calOneFace(const Eigen::Matrix<double, 3, 3> _threePoints)
		{
			Eigen::Matrix<double, 4, 1> result;
			result(0) = (_threePoints(1, 1) - _threePoints(1, 0)) * (_threePoints(2, 2) - _threePoints(2, 0)) - (_threePoints(2, 1) - _threePoints(2, 0)) * (_threePoints(1, 2) - _threePoints(1, 0));
			result(1) = (_threePoints(2, 1) - _threePoints(2, 0)) * (_threePoints(0, 2) - _threePoints(0, 0)) - (_threePoints(0, 1) - _threePoints(0, 0)) * (_threePoints(2, 2) - _threePoints(2, 0));
			result(2) = (_threePoints(0, 1) - _threePoints(0, 0)) * (_threePoints(1, 2) - _threePoints(1, 0)) - (_threePoints(1, 1) - _threePoints(1, 0)) * (_threePoints(0, 2) - _threePoints(0, 0));
			result(3) = -(result(0) * _threePoints(0, 0) + result(1) * _threePoints(1, 0) + result(2) * _threePoints(2, 0));
			return result;
		}

	public:
		SlopeEst(){
			points.resize(3, 4);
			estFaces.resize(4, 4);
			estFacesNormal.resize(3, 4);
		}

		/* 对类型的参数进行初始化 */
		void init(const Eigen::MatrixXd& _initPoints, const Eigen::Vector4i& _usedPoints)
		{
			this->points = _initPoints;
			this->usedPoints = _usedPoints;
			this->estRun();
		}

		/* 根据触地足有选择地进行更新 */
		void updatePoints(const Eigen::MatrixXd& _points, const Eigen::Vector4i& _contact)
		{
			for (int i = 0; i < 4; i++)
			{
				if (_contact(i) == 1)
				{
					this->points.col(i) = _points.col(i);
				}
			}
		}

		/* 计算可预测的平面，并得到各平面法向量的和向量 */
		const Eigen::Vector3d& estRun()
		{
			this->_calFaces();
			Eigen::Matrix<double, 3, 1> temp_normal;
			temp_normal.setZero();
			for (int i = 0; i < 4; i++)
			{
				/*double vel = pow(this->estFacesNormal(0, i), 2) + pow(this->estFacesNormal(1, i), 2) + pow(this->estFacesNormal(2, i), 2);
				if (vel > 0)
				{*/
					/*temp_normal(0) += this->estFacesNormal(0, i) / (4 * sqrt(vel));
					temp_normal(1) += this->estFacesNormal(1, i) / (4 * sqrt(vel));
					temp_normal(2) += this->estFacesNormal(2, i) / (4 * sqrt(vel));*/
					temp_normal += this->estFacesNormal.col(i).normalized();
				//}
			}
			this->estNormal = temp_normal.normalized();
			return this->estNormal;
		}

		/* 获取当前预测平面法向量 */
		const Eigen::Vector3d& getEstNormal()
		{
			return this->estNormal;
		}

		/* 通过估计的法向量计算估计平面的旋转矩阵 */
		Eigen::Matrix3d getSlopeRotation()
		{
			Vector3d refNormal(0, 0, 1);
			//Vector3d cross = refNormal.cross(this->estNormal);
			////std::cout << "axis: \n" << cross << std::endl;
			//double costheta = refNormal.dot(this->estNormal);
			//double sintheta = cross.norm();
			//cross.stableNormalize();
			///*std::cout << "axis: \n" << cross << std::endl;*/
			//Matrix3d K;
			//K << 0, -cross(2), cross(1), cross(2), 0, -cross(0), -cross(1), cross(0), 0;
			////std::cout << "K: \n" << K << std::endl;
			//Matrix3d R = Matrix3d::Identity() + sintheta * K + (1 - costheta) * K * K;
			//return R;
			return this->computeRotationMatrix(refNormal, this->estNormal);
		}

		Eigen::Matrix3d computeRotationMatrix(const Eigen::Vector3d& u, const Eigen::Vector3d& v) {
			// 归一化输入向量
			Eigen::Vector3d u_norm = u.normalized();
			Eigen::Vector3d v_norm = v.normalized();

			// 计算旋转轴
			Eigen::Vector3d k = u_norm.cross(v_norm);
			double k_norm = k.norm();

			// 处理共线情况（旋转角为 0 或 180°）
			if (k_norm < 1e-3) {
				// 如果 u 和 v 方向相同，返回单位矩阵
				if (u_norm.dot(v_norm) > 0) {
					return Eigen::Matrix3d::Identity();
				}
				// 如果 u 和 v 方向相反，返回绕任意垂直轴的 180° 旋转
				else {
					// 找到一个与 u 垂直的向量
					Eigen::Vector3d perpendicular = u_norm.unitOrthogonal();
					Eigen::AngleAxisd rotation(M_PI, perpendicular);
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
	};
}