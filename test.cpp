#include "Body.h"
#include "BodyCtrl.h"
#include "GaitCtrl.h"
#include "Leg.h"
#include "LegCtrl.h"

using namespace Quadruped;

int main()
{
    // self controll classes
    Leg lf_leg_obj(1, 1, 1, 1);
    LegCtrl lf_leg_ctrl(&lf_leg_obj, 2);
    lf_leg_ctrl.setEndPositionTar(Eigen::Vector3d(0, 0.0838, -0.27));

    Leg rf_leg_obj(1, 1, 1, 1);
    LegCtrl rf_leg_ctrl(&rf_leg_obj, 2);
    rf_leg_ctrl.setEndPositionTar(Eigen::Vector3d(0, -0.0838, -0.27));

    Leg lb_leg_obj(1, 1, 1, 1);
    LegCtrl lb_leg_ctrl(&lb_leg_obj, 2);
    lb_leg_ctrl.setEndPositionTar(Eigen::Vector3d(0, 0.0838, -0.27));

    Leg rb_leg_obj(1, 1, 1, 1);
    LegCtrl rb_leg_ctrl(&rb_leg_obj, 2);
    rb_leg_ctrl.setEndPositionTar(Eigen::Vector3d(0, -0.0838, -0.27));

    Leg* legsObj[4] = { &lf_leg_obj, &rf_leg_obj, &lb_leg_obj, &rb_leg_obj };
    Body qp_body(legsObj, static_cast<double>(2) * 0.001f);

    LegCtrl* legsCtrl[4] = { &lf_leg_ctrl,&rf_leg_ctrl,&lb_leg_ctrl,&rb_leg_ctrl };

    BodyCtrl qp_ctrl(&qp_body, legsCtrl, 2);
    bool use_mpc = false;

    Vector4d phaseResult;
    Vector4i contactResult;
    phaseResult.setZero();
    contactResult.setZero();
    GaitCtrl gaitCtrl(&qp_ctrl,legsCtrl,2,&phaseResult,&contactResult);
    
    return 0;
}