#pragma once

#include <Eigen/Dense>


namespace Quadruped {
	class SlopeEst {
	public:
		Eigen::MatrixXd points;// ���ڹ���ƽ��ĵ�
		Eigen::MatrixXd estFaces;// ����õ���ԭʼƽ��
		Eigen::MatrixXd estFacesNormal;// ����õ���ԭʼƽ�淨����
		Eigen::Vector3d estNormal;// ����ƽ��õ��Ĺ�һ��������
		Eigen::Vector4i usedPoints;// ����ʱʹ�õĵ�ʹ������

		/* ��������������ƽ�� */
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
				// ������������ƽ�����������ֱ�Ӷ���Ϊ������ֱ���ϵ�ƽ��
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
				// ��֤��������������
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
				// ע��Ϊ��ʹ��õ�ƽ�淨������ͬ���˴��Ķ��������ת����Ӧ��һ�£��˴�Ϊ��ʱ��
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
					// ��֤����������Ϊ����
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
			// ����õ��ĸ�ƽ��֮��ѷ��������
			for (int i = 0; i < 4; i++)
			{
				this->estFacesNormal = this->estFaces.block(0, 0, 3, 4);
			}
		}

		/* ͨ���������������ƽ��Ĳ��� */
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

		/* �����͵Ĳ������г�ʼ�� */
		void init(const Eigen::MatrixXd& _initPoints, const Eigen::Vector4i& _usedPoints)
		{
			this->points = _initPoints;
			this->usedPoints = _usedPoints;
			this->estRun();
		}

		/* ���ݴ�������ѡ��ؽ��и��� */
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

		/* �����Ԥ���ƽ�棬���õ���ƽ�淨�����ĺ����� */
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

		/* ��ȡ��ǰԤ��ƽ�淨���� */
		const Eigen::Vector3d& getEstNormal()
		{
			return this->estNormal;
		}

		/* ͨ�����Ƶķ������������ƽ�����ת���� */
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
			// ��һ����������
			Eigen::Vector3d u_norm = u.normalized();
			Eigen::Vector3d v_norm = v.normalized();

			// ������ת��
			Eigen::Vector3d k = u_norm.cross(v_norm);
			double k_norm = k.norm();

			// �������������ת��Ϊ 0 �� 180�㣩
			if (k_norm < 1e-3) {
				// ��� u �� v ������ͬ�����ص�λ����
				if (u_norm.dot(v_norm) > 0) {
					return Eigen::Matrix3d::Identity();
				}
				// ��� u �� v �����෴�����������ⴹֱ��� 180�� ��ת
				else {
					// �ҵ�һ���� u ��ֱ������
					Eigen::Vector3d perpendicular = u_norm.unitOrthogonal();
					Eigen::AngleAxisd rotation(M_PI, perpendicular);
					return rotation.toRotationMatrix();
				}
			}

			// ��һ����ת��
			k.normalize();

			// ������ת��
			double cos_theta = u_norm.dot(v_norm);
			double theta = acos(cos_theta);

			// ������Ԫ��
			Eigen::Quaterniond q;
			q = Eigen::AngleAxisd(theta, k);

			// ת��Ϊ��ת����
			return q.toRotationMatrix();
		}
	};
}