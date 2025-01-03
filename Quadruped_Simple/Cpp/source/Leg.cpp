#include "Leg.h"
#include <math.h>
#include <iostream>
#include "mathTool.h"

/* 初始化腿参数 */
void InitLeg(LegS* _leg, double* _links, double _ratio, int _dir)
{
#if DOF==3
	_leg->legPara.L[0] = _ratio * _links[0];
	_leg->legPara.L[1] = _links[1];
	_leg->legPara.L[2] = _links[2];
#elif DOF==2
	_leg->legPara.L[0] = _links[0];
	_leg->legPara.L[1] = _links[1];
#endif

	_leg->legPara.ratio = _ratio;
	_leg->legPara.dir = _dir;
}

/* 计算雅可比矩阵 */
void LegJacobiCal(LegS* _leg, JointS* _joint)
{
#if DOF==3
	// 正向雅可比（正运动学求微分）
	double s1 = sin(_joint->Angle[0]);
	double c1 = cos(_joint->Angle[0]);
	double s2 = sin(_joint->Angle[1]);
	double c2 = cos(_joint->Angle[1]);
	double s23 = sin(_joint->Angle[1] + _joint->Angle[2]);
	double c23 = cos(_joint->Angle[1] + _joint->Angle[2]);
	_leg->jacobi[0][0] = 0;
	_leg->jacobi[0][1] = -_leg->legPara.L[1] * c2 - _leg->legPara.L[2] * c23;
	_leg->jacobi[0][2] = -_leg->legPara.L[2] * c23;
	_leg->jacobi[1][0] = -_leg->legPara.L[0] * s1 + _leg->legPara.L[1] * c1 * c2 + _leg->legPara.L[2] * c1 * c23;
	_leg->jacobi[1][1] = -_leg->legPara.L[1] * s1 * s2 - _leg->legPara.L[2] * s1 * s23;
	_leg->jacobi[1][2] = -_leg->legPara.L[2] * s1 * s23;
	_leg->jacobi[2][0] = _leg->legPara.L[0] * c1 + _leg->legPara.L[1] * s1 * c2 + _leg->legPara.L[2] * s1 * c23;
	_leg->jacobi[2][1] = _leg->legPara.L[1] * c1 * s2 + _leg->legPara.L[2] * c1 * s23;
	_leg->jacobi[2][2] = _leg->legPara.L[2] * c1 * s23;
#elif DOF==2
	// 正向雅可比（正运动学求微分）
	double s1 = sin(_joint->Angle[0]);
	double c1 = cos(_joint->Angle[0]);
	double s2 = sin(_joint->Angle[1]);
	double c2 = cos(_joint->Angle[1]);
	double s12 = sin(_joint->Angle[0] + _joint->Angle[1]);
	double c12 = cos(_joint->Angle[0] + _joint->Angle[1]);
	_leg->jacobi[0][0] = -_leg->legPara.L[0] * c1 - _leg->legPara.L[1] * c12;
	_leg->jacobi[0][1] = -_leg->legPara.L[1] * c12;
	_leg->jacobi[2][0] = _leg->legPara.L[0] * s1 + _leg->legPara.L[1] * s12;
	_leg->jacobi[2][1] = _leg->legPara.L[1] * s12;
	m3d_transpose2(_leg->jacobiT, _leg->jacobi);
	
#endif
}

/* 计算逆雅可比矩阵 */
void LegInvJacobiCal(LegS* _leg, JointS* _joint)
{
#if DOF==3
	double l1 = _leg->legPara.L[0];
	double l2 = _leg->legPara.L[1];
	double l3 = _leg->legPara.L[2];
	double s1 = sin(_joint->Angle[0]);
	double c1 = cos(_joint->Angle[0]);
	double s2 = sin(_joint->Angle[1]);
	double c2 = cos(_joint->Angle[1]);
	double s23 = sin(_joint->Angle[1] + _joint->Angle[2]);
	double c23 = cos(_joint->Angle[1] + _joint->Angle[2]);
	double k1 = l2 * l3 * c23 * c1 * c1 * s2 - l2 * l3 * s23 * c1 * c1 * c2 + l2 * l3 * c23 * s1 * s1 * s2 - l2 * l3 * s23 * c2 * s1 * s1;
	double k2 = l3 * c23 + l2 * c2;
	double k3 = l1 * c1 + l3 * c23 * s1 + l2 * c2 * s1;
	double k4 = l3 * c23 * c1 - l1 * s1 + l2 * c1 * c2;
	double k5 = -s2 * l2 * l2 * c23 * c2 + s23 * l2 * l2 * c2 * c2
		- l3 * s2 * l2 * c23 * c23 + l3 * s23 * l2 * c23 * c2;
	_leg->jacobiI[0][0] = 0;
	_leg->jacobiI[0][1] = c1 / k2;
	_leg->jacobiI[0][2] = s1 / k2;
	_leg->jacobiI[1][0] = s23 / (l2 * c23 * s2 - l2 * s23 * c2);
	_leg->jacobiI[1][1] = c23 * k3 / k5;
	_leg->jacobiI[1][2] = c23 * k4 / k5;
	_leg->jacobiI[2][0] = (l3 * s23 + l2 * s2) / (l2 * l3 * c23 * s2 - l2 * l3 * s23 * c2);
	_leg->jacobiI[2][1] = k3 / k1;
	_leg->jacobiI[2][2] = k4 / k1;
#elif DOF==2
	double l1 = _leg->legPara.L[0];
	double l2 = _leg->legPara.L[1];
	double s1 = sin(_joint->Angle[0]);
	double c1 = cos(_joint->Angle[0]);
	double s2 = sin(_joint->Angle[1]);
	double c2 = cos(_joint->Angle[1]);
	double s12 = sin(_joint->Angle[0] + _joint->Angle[1]);
	double c12 = cos(_joint->Angle[0] + _joint->Angle[1]);
	double k1 = l1 * l2 * c12 * s1 - l1 * l2 * s12 * c1;
	double k2 = l1 * c12 * s1 - l1 * s12 * c1;
	_leg->jacobiI[0][0] = s12 / k2;
	_leg->jacobiI[0][2] = c12 / k2;
	_leg->jacobiI[2][0] = (l2 * s12 + l1 * s1) / k1;
	_leg->jacobiI[2][2] = (l2 * c12 + l1 * c1) / k1;
#endif
}

/* 正运动学计算 */
void LegFkCal(EndS* _end, JointS* _joint, LegParamS* _para)
{
#if DOF==3
	_end->Position[0] = -_para->L[1] * sin(_joint->Angle[1]) - _para->L[2] * sin(_joint->Angle[1] + _joint->Angle[2]);
	_end->Position[1] = _para->L[0] * cos(_joint->Angle[0]) + _para->L[1] * sin(_joint->Angle[0]) * cos(_joint->Angle[1]) + _para->L[2] * sin(_joint->Angle[0]) * cos(_joint->Angle[1] + _joint->Angle[2]);
	_end->Position[2] = _para->L[0] * sin(_joint->Angle[0]) - _para->L[1] * cos(_joint->Angle[0]) * cos(_joint->Angle[1]) - _para->L[2] * cos(_joint->Angle[0]) * cos(_joint->Angle[1] + _joint->Angle[2]);
#elif DOF==2
	_end->Position[0] = -_para->L[0] * sin(_joint->Angle[0]) - _para->L[1] * sin(_joint->Angle[0] + _joint->Angle[1]);
	_end->Position[2] = -_para->L[0] * cos(_joint->Angle[0]) - _para->L[1] * cos(_joint->Angle[0] + _joint->Angle[1]);
#endif
}

/* 逆运动学计算 */
void LegIkCal(JointS* _joint, EndS* _end, LegParamS* _para)
{
#if DOF==3
	double x = _end->Position[0];
	double y = _end->Position[1];
	double z = _end->Position[2];

	double L = sqrt(y * y + z * z - _para->L[0] * _para->L[0]);
	double theta1 = atan2((_para->L[0] * z + L * y), -L * z + _para->L[0] * y);
	L = sqrt(x * x + y * y + z * z - _para->L[0] * _para->L[0]);
	double theta3 = -3.1415926 + acos((_para->L[1] * _para->L[1] + _para->L[2] * _para->L[2] - L * L) / 2 / _para->L[1] / _para->L[2]);

	double a1 = y * sin(theta1) - z * cos(theta1);
	double m1 = -_para->L[2] * sin(theta3);
	double m2 = -_para->L[2] * cos(theta3) - _para->L[1];
	double theta2 = atan2(a1 * m1 + x * m2, x * m1 - a1 * m2);

	_joint->Angle[0] = theta1;
	_joint->Angle[1] = theta2;
	_joint->Angle[2] = theta3;
#elif DOF==2
	double x = _end->Position[0];
	double z = _end->Position[2];

	double L = sqrt(x * x + z * z);
	double theta = atan(x/z);
	//std::cout << "theta: " << theta << std::endl;
	double theta1 = 0;
	double theta2 = 0;
	if (_para->dir == -1)
	{
		theta1 = theta + acos((L * L + _para->L[0] * _para->L[0] - _para->L[1] * _para->L[1]) / (2 * _para->L[0] * L));
		theta2 = acos((-L * L + _para->L[0] * _para->L[0] + _para->L[1] * _para->L[1]) / (2 * _para->L[0] * _para->L[1])) - 3.1415926;
	}
	else if (_para->dir == 1)
	{
		theta1 = theta - acos((L * L + _para->L[0] * _para->L[0] - _para->L[1] * _para->L[1]) / (2 * _para->L[0] * L));
		theta2 = 3.1415926 - acos((-L * L + _para->L[0] * _para->L[0] + _para->L[1] * _para->L[1]) / (2 * _para->L[0] * _para->L[1]));
	}
	_joint->Angle[0] = theta1;
	_joint->Angle[1] = theta2;
#endif
}

