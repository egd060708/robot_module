#include "Quadruped/C/include/Body.h"
#include "Quadruped/C/include/GaitCtrl.h"
#include "Quadruped/C/include/Leg.h"
#include "Quadruped/C/include/LegCtrl.h"
#include <stdio.h>
#include <stdbool.h>
#include "mathPrint.h"

int main()
{
    /* 自定义类型 */
    LegS legsObj[4];
    Leg_Ctrl_Param legsCtrlParam[4];
    Body body;
    InitLeg(&legsObj[LF], 0.0838, 0.2, 0.2, 1);
    InitLeg(&legsObj[RF], 0.0838, 0.2, 0.2, -1);
    InitLeg(&legsObj[LB], 0.0838, 0.2, 0.2, 1);
    InitLeg(&legsObj[RB], 0.0838, 0.2, 0.2, -1);
    double kp_p[3] = { 10,10,-5 };
    //double kd_p[3] = { 0.5,0.5,-0.25 };
    double kd_p[3] = { 0.,0.,-0. };
    double kp_t[3] = { 1,1,1 };
    double init_p[4][3] = {
        {0, 0.0838, -0.27},
        {0, -0.0838, -0.27},
        {0, 0.0838, -0.27},
        {0, -0.0838, -0.27}
    };
    for (int i = 0; i < 4; i++)
    {
        InitLegCtrl(&legsCtrlParam[i], kp_p, kd_p, kp_t, (double)2 * 0.001);
        StartVirtualVMC(init_p[i], &legsObj[i]);
    }
    LegS* lobj[4] = { &legsObj[LF],&legsObj[RF],&legsObj[LB],&legsObj[RB] };
    InitBody(&body, lobj, (double)2 * 0.001);
    double leg2Body[3] = { 0.1805, 0.047, 0 };
    SetBodyParam(&body, leg2Body);
    double init_Ep[4][3] = {
        {0.1805, 0.1308, 0.},
        {0.1805, -0.1308, 0.},
        {-0.1805, 0.1308, 0.},
        {-0.1805, -0.1308, 0.},
    };
    UpdateFootPoint(&body.targetWorldState, init_Ep);
    UpdateFootPointBase(&body.targetWorldState, init_Ep);
    double angle_t[3] = { 0,0,0 };
    double p_t[3] = { 0,0,0.27 };
    UpdateBodyPos(&body.targetWorldState, angle_t, p_t);
    GaitS gaitData;
    double bias[4] = { 0.5, 0, 0, 0.5 };
    InitGaitCtrl(&gaitData, &body, 0.4, 0.5, bias);
    double cmd_slope[3] = { 0.01,0.01,0.01 };
    bool is_gait = false;
    printVd("cmd",cmd_slope,3);
    system("pause");
    return 0;
}