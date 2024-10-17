#include "water_gait1.h"

#define PI 3.1415926

// forelimb正
//double forelimb[ROWS][COLS] = {{-2.609164, 0.156027}, {-2.565879, 0.167102}, {-2.521802, 0.176457}, \
//        {-2.477029, 0.184087}, {-2.431657, 0.189988}, {-2.385776, 0.194160}, \
//        {-2.339479, 0.196605}, {-2.292850, 0.197328}, {-2.245976, 0.196339}, \
//        {-2.198937, 0.193647}, {-2.151810, 0.189266}, {-2.104669, 0.183211}, \
//        {-2.057585, 0.175501}, {-2.010625, 0.166155}, {-1.963852, 0.155197}, \
//        {-1.917324, 0.142652}, {-1.871098, 0.128544}, {-1.825224, 0.112905}, \
//        {-1.779751, 0.095762}, {-1.734723, 0.077150}, {-1.690181, 0.057102}, \
//        {-1.646161, 0.035654}, {-1.602698, 0.012842}, {-1.571117, 0.011294}, \
//        {-1.554279, 0.036715}, {-1.539324, 0.063379}, {-1.526232, 0.091243}, \
//        {-1.514980, 0.120263}, {-1.505541, 0.150394}, {-1.497884, 0.181589}, \
//        {-1.491977, 0.213800}, {-1.487784, 0.246979}, {-1.485267, 0.281075}, \
//        {-1.484388, 0.316039}, {-1.485105, 0.351817}, {-1.487377, 0.388358}, \
//        {-1.491160, 0.425606}, {-1.496410, 0.463509}, {-1.503083, 0.502009}, \
//        {-1.511135, 0.541051}, {-1.520521, 0.580577}, {-1.531198, 0.620529}, \
//        {-1.543120, 0.660849}, {-1.556246, 0.701477}, {-1.570533, 0.742354}, \
//        {-1.585940, 0.783417}, {-1.602428, 0.824607}, {-1.619957, 0.865861}, \
//        {-1.638490, 0.907118}, {-1.657993, 0.948314}, {-1.678430, 0.989386}, \
//        {-1.699770, 1.030272}, {-1.721981, 1.070909}, {-1.745035, 1.111232}, \
//        {-1.768904, 1.151178}, {-1.793563, 1.190684}, {-1.818987, 1.229686}, \
//        {-1.845154, 1.268122}, {-1.872042, 1.305928}, {-1.899633, 1.343043}, \
//        {-1.927907, 1.379405}, {-1.956847, 1.414953}, {-1.986437, 1.449627}, \
//        {-2.016661, 1.483368}, {-2.047505, 1.516117}, {-2.078955, 1.547818}, \
//        {-2.110996, 1.578415}, {-2.143615, 1.607854}, {-2.176799, 1.636082}, \
//        {-2.210533, 1.663049}, {-2.244803, 1.688705}, {-2.279594, 1.713003}, \
//        {-2.314891, 1.735898}, {-2.350677, 1.757349}, {-2.386934, 1.777313}, \
//        {-2.423643, 1.795753}, {-2.460784, 1.812634}, {-2.498334, 1.827923}, \
//        {-2.536269, 1.841591}, {-2.574564, 1.853609}, {-2.613188, 1.863955}, \
//        {-2.652114, 1.872607}, {-2.691306, 1.879548}, {-2.730731, 1.884762}, \
//        {-2.770351, 1.888239}, {-2.810124, 1.889970}, {-2.850008, 1.889951}, \
//        {-2.889958, 1.888181}, {-2.929926, 1.884661}, {-2.969859, 1.879398}, \
//        {-3.009706, 1.872400}, {-3.049409, 1.863679}, {-3.088910, 1.853253}, \
//        {-3.128149, 1.841138}, {3.116123, 1.827359}, {3.077602, 1.811941}, \
//        {3.039539, 1.794912}, {3.002003, 1.776304}, {2.965067, 1.756153}, \
//        {2.928802, 1.734497}, {2.893282, 1.711376}, {2.858582, 1.686834}, \
//        {2.824777, 1.660917}, {2.791943, 1.633673}, {2.760157, 1.605154}, \
//        {2.729496, 1.575412}, {2.700035, 1.544503}, {2.671850, 1.512483}, \
//        {2.645016, 1.479412}, {2.619607, 1.445350}, {2.595695, 1.410358}, \
//        {2.573351, 1.374500}, {2.552642, 1.337839}, {2.533637, 1.300442}, \
//        {2.516397, 1.262373}, {2.500985, 1.223699}, {2.487457, 1.184487}, \
//        {2.475868, 1.144804}, {2.466267, 1.104718}, {2.458702, 1.064295}, \
//        {2.453215, 1.023603}, {2.449844, 0.982708}, {2.448621, 0.941676}, \
//        {2.449574, 0.900574}, {2.452728, 0.859466}, {2.458100, 0.818417}, \
//        {2.465702, 0.777489}, {2.475542, 0.736744}, {2.487622, 0.696245}, \
//        {2.501937, 0.656050}, {2.518477, 0.616218}, {2.537227, 0.576807}, \
//        {2.558165, 0.537872}, {2.581264, 0.499466}, {2.606492, 0.461644}, \
//        {2.633809, 0.424457}, {2.663172, 0.387953}, {2.694530, 0.352181}, \
//        {2.727828, 0.317188}, {2.763005, 0.283019}, {2.799996, 0.249716}, \
//        {2.838729, 0.217321}, {2.879130, 0.185875}, {2.921118, 0.155417}, \
//        {2.964609, 0.125982}, {3.009514, 0.097608}, {3.055741, 0.070327}, \
//        {3.103194, 0.044172}, {-3.131412, 0.019175}, {3.093238, 0.135772}, \
//        {-3.134655, 0.104372}, {-3.078635, 0.074316}, {-3.022025, 0.045653}, \
//        {-2.964963, 0.018428}, {-2.914903, 0.007316}, {-2.881573, 0.031542}, \
//        {-2.846645, 0.054213}, {-2.810218, 0.075297}, {-2.772394, 0.094766}, \
//        {-2.733275, 0.112594}, {-2.692963, 0.128759}, {-2.651559, 0.143242}};

// double forelimb_cont[ROWS][COLS] = {{-2.498334, 1.827923}, {-2.536269, 1.841591}, {-2.574564, 1.853609}, \
//         {-2.613188, 1.863955}, {-2.652114, 1.872607}, {-2.691306, 1.879548}, \
//         {-2.730731, 1.884762}, {-2.770351, 1.888239}, {-2.810124, 1.889970}, \
//         {-2.850008, 1.889951}, {-2.889958, 1.888181}, {-2.929926, 1.884661}, \
//         {-2.969859, 1.879398}, {-3.009706, 1.872400}, {-3.049409, 1.863679}, \
//         {-3.088910, 1.853253}, {-3.128149, 1.841138}, {3.116123, 1.827359}, \
//         {3.077602, 1.811941}, {3.039539, 1.794912}, {3.002003, 1.776304}, \
//         {2.965067, 1.756153}, {2.928802, 1.734497}, {2.893282, 1.711376}, \
//         {2.858582, 1.686834}, {2.824777, 1.660917}, {2.791943, 1.633673}, \
//         {2.760157, 1.605154}, {2.729496, 1.575412}, {2.700035, 1.544503}, \
//         {2.671850, 1.512483}, {2.645016, 1.479412}, {2.619607, 1.445350}, \
//         {2.595695, 1.410358}, {2.573351, 1.374500}, {2.552642, 1.337839}, \
//         {2.533637, 1.300442}, {2.516397, 1.262373}, {2.500985, 1.223699}, \
//         {2.487457, 1.184487}, {2.475868, 1.144804}, {2.466267, 1.104718}, \
//         {2.458702, 1.064295}, {2.453215, 1.023603}, {2.449844, 0.982708}, \
//         {2.448621, 0.941676}, {2.449574, 0.900574}, {2.452728, 0.859466}, \
//         {2.458100, 0.818417}, {2.465702, 0.777489}, {2.475542, 0.736744}, \
//         {2.487622, 0.696245}, {2.501937, 0.656050}, {2.518477, 0.616218}, \
//         {2.537227, 0.576807}, {2.558165, 0.537872}, {2.581264, 0.499466}, \
//         {2.606492, 0.461644}, {2.633809, 0.424457}, {2.663172, 0.387953}, \
//         {2.694530, 0.352181}, {2.727828, 0.317188}, {2.763005, 0.283019}, \
//         {2.799996, 0.249716}, {2.838730, 0.217321}, {2.879130, 0.185875}, \
//         {2.921118, 0.155417}, {2.964609, 0.125982}, {3.009514, 0.097608}, \
//         {3.055741, 0.070327}, {3.103194, 0.044172}, {-3.131412, 0.019175}, \
//         {3.093238, 0.135772}, {-3.134655, 0.104372}, {-3.078635, 0.074316}, \
//         {-3.022024, 0.045653}, {-2.964963, 0.018428}, {-2.914903, 0.007316}, \
//         {-2.881573, 0.031542}, {-2.844565, 0.050054}, {-2.810218, 0.075297}, \
//         {-2.772181, 0.094339}, {-2.733275, 0.112594}, {-2.692963, 0.128759}, \
//         {-2.651559, 0.143242}, {-2.609164, 0.156027}, {-2.565879, 0.167102}, \
//         {-2.521802, 0.176457}, {-2.477029, 0.184087}, {-2.431657, 0.189988}, \
//         {-2.385776, 0.194160}, {-2.339479, 0.196605}, {-2.292850, 0.197328}, \
//         {-2.245976, 0.196339}, {-2.198937, 0.193647}, {-2.151810, 0.189266}, \
//         {-2.104669, 0.183211}, {-2.057585, 0.175501}, {-2.010625, 0.166155}, \
//         {-1.963852, 0.155197}, {-1.917324, 0.142652}, {-1.871098, 0.128544}, \
//         {-1.825224, 0.112905}, {-1.779751, 0.095762}, {-1.734723, 0.077150}, \
//         {-1.690181, 0.057102}, {-1.646161, 0.035654}, {-1.602698, 0.012842}, \
//         {-1.571117, 0.011294}, {-1.554279, 0.036715}, {-1.539324, 0.063379}, \
//         {-1.526232, 0.091243}, {-1.514980, 0.120263}, {-1.505541, 0.150394}, \
//         {-1.497884, 0.181589}, {-1.491977, 0.213800}, {-1.487784, 0.246979}, \
//         {-1.485267, 0.281075}, {-1.484388, 0.316039}, {-1.485105, 0.351817}, \
//         {-1.487377, 0.388358}, {-1.491160, 0.425606}, {-1.496410, 0.463509}, \
//         {-1.503083, 0.502009}, {-1.511135, 0.541051}, {-1.520521, 0.580577}, \
//         {-1.531198, 0.620529}, {-1.543120, 0.660849}, {-1.556246, 0.701477}, \
//         {-1.570533, 0.742354}, {-1.585940, 0.783417}, {-1.602428, 0.824607}, \
//         {-1.619957, 0.865861}, {-1.638490, 0.907118}, {-1.657993, 0.948314}, \
//         {-1.678430, 0.989386}, {-1.699770, 1.030272}, {-1.721981, 1.070909}, \
//         {-1.745035, 1.111232}, {-1.768904, 1.151178}, {-1.793563, 1.190684}, \
//         {-1.818987, 1.229686}, {-1.845154, 1.268122}, {-1.872042, 1.305928}, \
//         {-1.899633, 1.343043}, {-1.927907, 1.379405}, {-1.956847, 1.414953}, \
//         {-1.986437, 1.449627}, {-2.016661, 1.483368}, {-2.047505, 1.516117}, \
//         {-2.078955, 1.547818}, {-2.110996, 1.578415}, {-2.143615, 1.607854}, \
//         {-2.176799, 1.636082}, {-2.210533, 1.663049}, {-2.244803, 1.688705}, \
//         {-2.279594, 1.713003}, {-2.314891, 1.735898}, {-2.350677, 1.757349}, \
//         {-2.386934, 1.777313}, {-2.423643, 1.795753}, {-2.460784, 1.812634}};

// // hindlimb反
// double hindlimb[ROWS][COLS] = {{-0.865582, -0.250746}, {-0.839609, -0.233767}, {-0.813066, -0.219054}, \
//         {-0.785986, -0.206649}, {-0.758404, -0.196587}, {-0.730357, -0.188892}, \
//         {-0.701890, -0.183580}, {-0.673047, -0.180655}, {-0.643880, -0.180115}, \
//         {-0.614441, -0.181945}, {-0.584784, -0.186124}, {-0.554969, -0.192619}, \
//         {-0.525056, -0.201389}, {-0.495106, -0.212383}, {-0.465185, -0.225543}, \
//         {-0.435356, -0.240801}, {-0.405687, -0.258083}, {-0.376244, -0.277305}, \
//         {-0.347094, -0.298378}, {-0.318305, -0.321206}, {-0.289942, -0.345686}, \
//         {-0.262071, -0.371710}, {-0.234758, -0.399167}, {-0.208065, -0.427940}, \
//         {-0.182054, -0.457908}, {-0.156784, -0.488949}, {-0.132313, -0.520938}, \
//         {-0.108694, -0.553749}, {-0.085980, -0.587255}, {-0.064220, -0.621330}, \
//         {-0.043457, -0.655846}, {-0.023734, -0.690679}, {-0.005090, -0.725707}, \
//         {0.012443, -0.760810}, {0.028834, -0.795872}, {0.044055, -0.830780}, \
//         {0.058087, -0.865425}, {0.070909, -0.899706}, {0.082509, -0.933525}, \
//         {0.092878, -0.966790}, {0.102008, -0.999418}, {0.109899, -1.031329}, \
//         {0.116554, -1.062453}, {0.121978, -1.092727}, {0.126181, -1.122094}, \
//         {0.129178, -1.150506}, {0.130984, -1.177922}, {0.131621, -1.204310}, \
//         {0.131110, -1.229644}, {0.129478, -1.253908}, {0.126754, -1.277092}, \
//         {0.122968, -1.299193}, {0.118152, -1.320217}, {0.112340, -1.340175}, \
//         {0.105570, -1.359086}, {0.097877, -1.376974}, {0.089299, -1.393869}, \
//         {0.079876, -1.409808}, {0.069646, -1.424830}, {0.058648, -1.438980}, \
//         {0.046922, -1.452308}, {0.034506, -1.464865}, {0.021437, -1.476705}, \
//         {0.007754, -1.487884}, {-0.006508, -1.498461}, {-0.021314, -1.508494}, \
//         {-0.036630, -1.518042}, {-0.052424, -1.527163}, {-0.068666, -1.535916}, \
//         {-0.085328, -1.544354}, {-0.102382, -1.552533}, {-0.119805, -1.560503}, \
//         {-0.137574, -1.568311}, {-0.155667, -1.576002}, {-0.174067, -1.583616}, \
//         {-0.192756, -1.591187}, {-0.211721, -1.598747}, {-0.230948, -1.606319}, \
//         {-0.250427, -1.613923}, {-0.270148, -1.621574}, {-0.290104, -1.629277}, \
//         {-0.310289, -1.637035}, {-0.330697, -1.644843}, {-0.351325, -1.652690}, \
//         {-0.372170, -1.660558}, {-0.393229, -1.668424}, {-0.414502, -1.676257}, \
//         {-0.435986, -1.684024}, {-0.457681, -1.691681}, {-0.479584, -1.699184}, \
//         {-0.501694, -1.706479}, {-0.524009, -1.713512}, {-0.546525, -1.720221}, \
//         {-0.569237, -1.726541}, {-0.592140, -1.732406}, {-0.615226, -1.737745}, \
//         {-0.638487, -1.742484}, {-0.661912, -1.746549}, {-0.685487, -1.749864}, \
//         {-0.709197, -1.752354}, {-0.733025, -1.753941}, {-0.756952, -1.754551}, \
//         {-0.780953, -1.754109}, {-0.805004, -1.752544}, {-0.829077, -1.749785}, \
//         {-0.853142, -1.745767}, {-0.877164, -1.740427}, {-0.901107, -1.733707}, \
//         {-0.924932, -1.725554}, {-0.948598, -1.715920}, {-0.972060, -1.704763}, \
//         {-0.995272, -1.692048}, {-1.018185, -1.677746}, {-1.040747, -1.661837}, \
//         {-1.062906, -1.644306}, {-1.084606, -1.625146}, {-1.105791, -1.604362}, \
//         {-1.126403, -1.581961}, {-1.146384, -1.557964}, {-1.165675, -1.532397}, \
//         {-1.184215, -1.505295}, {-1.201944, -1.476703}, {-1.218802, -1.446673}, \
//         {-1.234731, -1.415264}, {-1.249671, -1.382547}, {-1.263565, -1.348595}, \
//         {-1.276357, -1.313494}, {-1.287994, -1.277334}, {-1.298424, -1.240211}, \
//         {-1.307597, -1.202230}, {-1.315468, -1.163499}, {-1.321992, -1.124134}, \
//         {-1.327131, -1.084253}, {-1.330848, -1.043979}, {-1.333111, -1.003438}, \
//         {-1.333894, -0.962761}, {-1.333172, -0.922079}, {-1.330928, -0.881524}, \
//         {-1.327148, -0.841231}, {-1.321824, -0.801333}, {-1.314953, -0.761963}, \
//          {-1.306538, -0.723252}, {-1.296586, -0.685332}, {-1.285110, -0.648329}, \
//         {-1.272131, -0.612366}, {-1.257671, -0.577563}, {-1.241761, -0.544036}, \
//         {-1.224436, -0.511894}, {-1.205737, -0.481242}, {-1.185710, -0.452177}, \
//         {-1.164405, -0.424792}, {-1.141879, -0.399168}, {-1.118192, -0.375384}, \
//         {-1.093409, -0.353508}, {-1.067599, -0.333600}, {-1.040834, -0.315712}, \
//         {-1.013192, -0.299888}, {-0.984752, -0.286161}, {-0.955597, -0.274557}, \
//         {-0.925812, -0.265094}, {-0.895484, -0.257779}, {-0.864703, -0.252611}};

// double hindlimb_cont[ROWS][COLS] = {{-0.351325, -1.652690}, {-0.372170, -1.660558}, {-0.393229, -1.668424}, \
//         {-0.414502, -1.676257}, {-0.435986, -1.684024}, {-0.457681, -1.691681}, \
//         {-0.479584, -1.699184}, {-0.501694, -1.706479}, {-0.524009, -1.713512}, \
//         {-0.546525, -1.720221}, {-0.569237, -1.726542}, {-0.592140, -1.732406}, \
//         {-0.615226, -1.737745}, {-0.638487, -1.742484}, {-0.661912, -1.746550}, \
//         {-0.685487, -1.749864}, {-0.709197, -1.752354}, {-0.733025, -1.753941}, \
//         {-0.756952, -1.754551}, {-0.780953, -1.754109}, {-0.805004, -1.752544}, \
//         {-0.829077, -1.749785}, {-0.853142, -1.745767}, {-0.877164, -1.740427}, \
//         {-0.901107, -1.733707}, {-0.924932, -1.725554}, {-0.948598, -1.715920}, \
//         {-0.972060, -1.704763}, {-0.995272, -1.692048}, {-1.018185, -1.677746}, \
//         {-1.040747, -1.661837}, {-1.062906, -1.644306}, {-1.084606, -1.625146}, \
//         {-1.105791, -1.604362}, {-1.126403, -1.581961}, {-1.146384, -1.557964}, \
//         {-1.165675, -1.532397}, {-1.184215, -1.505295}, {-1.201944, -1.476703}, \
//         {-1.218802, -1.446673}, {-1.234731, -1.415264}, {-1.249671, -1.382547}, \
//         {-1.263565, -1.348595}, {-1.276357, -1.313494}, {-1.287994, -1.277334}, \
//         {-1.298424, -1.240211}, {-1.307597, -1.202230}, {-1.315468, -1.163499}, \
//         {-1.321992, -1.124134}, {-1.327131, -1.084253}, {-1.330848, -1.043979}, \
//         {-1.333111, -1.003438}, {-1.333894, -0.962761}, {-1.333172, -0.922079}, \
//         {-1.330928, -0.881524}, {-1.327148, -0.841231}, {-1.321824, -0.801333}, \
//         {-1.314953, -0.761963}, {-1.306538, -0.723252}, {-1.296586, -0.685332}, \
//         {-1.285110, -0.648329}, {-1.272131, -0.612366}, {-1.257671, -0.577563}, \
//         {-1.241761, -0.544036}, {-1.224436, -0.511894}, {-1.205737, -0.481242}, \
//         {-1.185710, -0.452177}, {-1.164405, -0.424792}, {-1.141879, -0.399168}, \
//         {-1.118192, -0.375384}, {-1.093409, -0.353508}, {-1.067599, -0.333600}, \
//         {-1.040834, -0.315712}, {-1.013192, -0.299888}, {-0.984752, -0.286161}, \
//         {-0.955597, -0.274557}, {-0.925812, -0.265094}, {-0.895484, -0.257779}, \
//         {-0.864703, -0.252611}, {-0.865582, -0.250746}, {-0.839609, -0.233767}, \
//         {-0.813066, -0.219054}, {-0.785986, -0.206649}, {-0.758404, -0.196587}, \
//         {-0.730357, -0.188892}, {-0.701890, -0.183580}, {-0.673047, -0.180655}, \
//         {-0.643880, -0.180115}, {-0.614441, -0.181945}, {-0.584784, -0.186124}, \
//         {-0.554969, -0.192619}, {-0.525056, -0.201389}, {-0.495106, -0.212383}, \
//         {-0.465185, -0.225543}, {-0.435356, -0.240801}, {-0.405687, -0.258083}, \
//         {-0.376244, -0.277305}, {-0.347094, -0.298378}, {-0.318305, -0.321206}, \
//         {-0.289942, -0.345686}, {-0.262071, -0.371710}, {-0.234758, -0.399167}, \
//         {-0.208065, -0.427940}, {-0.182054, -0.457908}, {-0.156784, -0.488949}, \
//         {-0.132313, -0.520938}, {-0.108694, -0.553749}, {-0.085980, -0.587255}, \
//         {-0.064220, -0.621330}, {-0.043457, -0.655846}, {-0.023734, -0.690679}, \
//         {-0.005090, -0.725707}, {0.012443, -0.760810}, {0.028834, -0.795872}, \
//         {0.044055, -0.830780}, {0.058087, -0.865425}, {0.070909, -0.899706}, \
//         {0.082509, -0.933525}, {0.092878, -0.966790}, {0.102008, -0.999418}, \
//         {0.109899, -1.031329}, {0.116554, -1.062453}, {0.121978, -1.092727}, \
//         {0.126181, -1.122094}, {0.129178, -1.150506}, {0.130984, -1.177922}, \
//         {0.131621, -1.204311}, {0.131110, -1.229644}, {0.129478, -1.253908}, \
//         {0.126754, -1.277092}, {0.122968, -1.299193}, {0.118152, -1.320217}, \
//         {0.112340, -1.340175}, {0.105570, -1.359086}, {0.097877, -1.376974}, \
//         {0.089299, -1.393869}, {0.079876, -1.409808}, {0.069646, -1.424830}, \
//         {0.058648, -1.438980}, {0.046922, -1.452308}, {0.034506, -1.464865}, \
//         {0.021437, -1.476705}, {0.007754, -1.487884}, {-0.006508, -1.498461}, \
//         {-0.021314, -1.508494}, {-0.036630, -1.518042}, {-0.052424, -1.527163}, \
//         {-0.068666, -1.535916}, {-0.085328, -1.544354}, {-0.102382, -1.552533}, \
//         {-0.119805, -1.560503}, {-0.137574, -1.568311}, {-0.155667, -1.576002}, \
//         {-0.174067, -1.583616}, {-0.192756, -1.591187}, {-0.211721, -1.598747}, \
//         {-0.230948, -1.606319}, {-0.250427, -1.613923}, {-0.270148, -1.621574}, \
//         {-0.290104, -1.629277}, {-0.310289, -1.637035}, {-0.330697, -1.644843}};

double forelimb[ROWS][COLS] = {{-2.609164, 0.156027}, {-2.565879, 0.167102}, {-2.521802, 0.176457}, \
       {-2.477029, 0.184087}, {-2.431657, 0.189988}, {-2.385776, 0.194160}, \
       {-2.339479, 0.196605}, {-2.292850, 0.197328}, {-2.245976, 0.196339}, \
       {-2.198937, 0.193647}, {-2.151810, 0.189266}, {-2.104669, 0.183211}, \
       {-2.057585, 0.175501}, {-2.010625, 0.166155}, {-1.963852, 0.155197}, \
       {-1.917324, 0.142652}, {-1.871098, 0.128544}, {-1.825224, 0.112905}, \
       {-1.779751, 0.095762}, {-1.734723, 0.077150}, {-1.690181, 0.057102}, \
       {-1.646161, 0.035654}, {-1.602698, 0.012842}, {-1.571117, 0.011294}, \
       {-1.554279, 0.036715}, {-1.539324, 0.063379}, {-1.526232, 0.091243}, \
       {-1.514980, 0.120263}, {-1.505541, 0.150394}, {-1.497884, 0.181589}, \
       {-1.491977, 0.213800}, {-1.487784, 0.246979}, {-1.485267, 0.281075}, \
       {-1.484388, 0.316039}, {-1.485105, 0.351817}, {-1.487377, 0.388358}, \
       {-1.491160, 0.425606}, {-1.496410, 0.463509}, {-1.503083, 0.502009}, \
       {-1.511135, 0.541051}, {-1.520521, 0.580577}, {-1.531198, 0.620529}, \
       {-1.543120, 0.660849}, {-1.556246, 0.701477}, {-1.570533, 0.742354}, \
       {-1.585940, 0.783417}, {-1.602428, 0.824607}, {-1.619957, 0.865861}, \
       {-1.638490, 0.907118}, {-1.657993, 0.948314}, {-1.678430, 0.989386}, \
       {-1.699770, 1.030272}, {-1.721981, 1.070909}, {-1.745035, 1.111232}, \
       {-1.768904, 1.151178}, {-1.793563, 1.190684}, {-1.818987, 1.229686}, \
       {-1.845154, 1.268122}, {-1.872042, 1.305928}, {-1.899633, 1.343043}, \
       {-1.927907, 1.379405}, {-1.956847, 1.414953}, {-1.986437, 1.449627}, \
       {-2.016661, 1.483368}, {-2.047505, 1.516117}, {-2.078955, 1.547818}, \
       {-2.110996, 1.578415}, {-2.143615, 1.607854}, {-2.176799, 1.636082}, \
       {-2.210533, 1.663049}, {-2.244803, 1.688705}, {-2.279594, 1.713003}, \
       {-2.314891, 1.735898}, {-2.350677, 1.757349}, {-2.386934, 1.777313}, \
       {-2.423643, 1.795753}, {-2.460784, 1.812634}, {-2.498334, 1.827923}, \
       {-2.536269, 1.841591}, {-2.574564, 1.853609}, {-2.613188, 1.863955}, \
       {-2.652114, 1.872607}, {-2.691306, 1.879548}, {-2.730731, 1.884762}, \
       {-2.770351, 1.888239}, {-2.810124, 1.889970}, {-2.850008, 1.889951}, \
       {-2.889958, 1.888181}, {-2.929926, 1.884661}, {-2.969859, 1.879398}, \
       {-3.009706, 1.872400}, {-3.049409, 1.863679}, {-3.088910, 1.853253}, \
       {-3.128149, 1.841138}, {3.116123-2*PI, 1.827359}, {3.077602-2*PI, 1.811941}, \
       {3.039539-2*PI, 1.794912}, {3.002003-2*PI, 1.776304}, {2.965067-2*PI, 1.756153}, \
       {2.928802-2*PI, 1.734497}, {2.893282-2*PI, 1.711376}, {2.858582-2*PI, 1.686834}, \
       {2.824777-2*PI, 1.660917}, {2.791943-2*PI, 1.633673}, {2.760157-2*PI, 1.605154}, \
       {2.729496-2*PI, 1.575412}, {2.700035-2*PI, 1.544503}, {2.671850-2*PI, 1.512483}, \
       {2.645016-2*PI, 1.479412}, {2.619607-2*PI, 1.445350}, {2.595695-2*PI, 1.410358}, \
       {2.573351-2*PI, 1.374500}, {2.552642-2*PI, 1.337839}, {2.533637-2*PI, 1.300442}, \
       {2.516397-2*PI, 1.262373}, {2.500985-2*PI, 1.223699}, {2.487457-2*PI, 1.184487}, \
       {2.475868-2*PI, 1.144804}, {2.466267-2*PI, 1.104718}, {2.458702-2*PI, 1.064295}, \
       {2.453215-2*PI, 1.023603}, {2.449844-2*PI, 0.982708}, {2.448621-2*PI, 0.941676}, \
       {2.449574-2*PI, 0.900574}, {2.452728-2*PI, 0.859466}, {2.458100-2*PI, 0.818417}, \
       {2.465702-2*PI, 0.777489}, {2.475542-2*PI, 0.736744}, {2.487622-2*PI, 0.696245}, \
       {2.501937-2*PI, 0.656050}, {2.518477-2*PI, 0.616218}, {2.537227-2*PI, 0.576807}, \
       {2.558165-2*PI, 0.537872}, {2.581264-2*PI, 0.499466}, {2.606492-2*PI, 0.461644}, \
       {2.633809-2*PI, 0.424457}, {2.663172-2*PI, 0.387953}, {2.694530-2*PI, 0.352181}, \
       {2.727828-2*PI, 0.317188}, {2.763005-2*PI, 0.283019}, {2.799996-2*PI, 0.249716}, \
       {2.838729-2*PI, 0.217321}, {2.879130-2*PI, 0.185875}, {2.921118-2*PI, 0.155417}, \
       {2.964609-2*PI, 0.125982}, {3.009514-2*PI, 0.097608}, {3.055741-2*PI, 0.070327}, \
       {3.103194-2*PI, 0.044172}, {-3.131412, 0.019175}, {3.093238 - 2 * PI, 0.135772}, \
       {-3.134655, 0.104372}, {-3.078635, 0.074316}, {-3.022025, 0.045653}, \
       {-2.964963, 0.018428}, {-2.914903, 0.007316}, {-2.881573, 0.031542}, \
       {-2.846645, 0.054213}, {-2.810218, 0.075297}, {-2.772394, 0.094766}, \
       {-2.733275, 0.112594}, {-2.692963, 0.128759}, {-2.651559, 0.143242}};

double forelimb_cont[ROWS][COLS] = {{-2.498334, 1.827923}, {-2.536269, 1.841591}, {-2.574564, 1.853609}, \
        {-2.613188, 1.863955}, {-2.652114, 1.872607}, {-2.691306, 1.879548}, \
        {-2.730731, 1.884762}, {-2.770351, 1.888239}, {-2.810124, 1.889970}, \
        {-2.850008, 1.889951}, {-2.889958, 1.888181}, {-2.929926, 1.884661}, \
        {-2.969859, 1.879398}, {-3.009706, 1.872400}, {-3.049409, 1.863679}, \
        {-3.088910, 1.853253}, {-3.128149, 1.841138}, {3.116123 - 2 * PI, 1.827359}, \
        {3.077602-2*PI, 1.811941}, {3.039539-2*PI, 1.794912}, {3.002003-2*PI, 1.776304}, \
        {2.965067-2*PI, 1.756153}, {2.928802-2*PI, 1.734497}, {2.893282-2*PI, 1.711376}, \
        {2.858582-2*PI, 1.686834}, {2.824777-2*PI, 1.660917}, {2.791943-2*PI, 1.633673}, \
        {2.760157-2*PI, 1.605154}, {2.729496-2*PI, 1.575412}, {2.700035-2*PI, 1.544503}, \
        {2.671850-2*PI, 1.512483}, {2.645016-2*PI, 1.479412}, {2.619607-2*PI, 1.445350}, \
        {2.595695-2*PI, 1.410358}, {2.573351-2*PI, 1.374500}, {2.552642-2*PI, 1.337839}, \
        {2.533637-2*PI, 1.300442}, {2.516397-2*PI, 1.262373}, {2.500985-2*PI, 1.223699}, \
        {2.487457-2*PI, 1.184487}, {2.475868-2*PI, 1.144804}, {2.466267-2*PI, 1.104718}, \
        {2.458702-2*PI, 1.064295}, {2.453215-2*PI, 1.023603}, {2.449844-2*PI, 0.982708}, \
        {2.448621-2*PI, 0.941676}, {2.449574-2*PI, 0.900574}, {2.452728-2*PI, 0.859466}, \
        {2.458100-2*PI, 0.818417}, {2.465702-2*PI, 0.777489}, {2.475542-2*PI, 0.736744}, \
        {2.487622-2*PI, 0.696245}, {2.501937-2*PI, 0.656050}, {2.518477-2*PI, 0.616218}, \
        {2.537227-2*PI, 0.576807}, {2.558165-2*PI, 0.537872}, {2.581264-2*PI, 0.499466}, \
        {2.606492-2*PI, 0.461644}, {2.633809-2*PI, 0.424457}, {2.663172-2*PI, 0.387953}, \
        {2.694530-2*PI, 0.352181}, {2.727828-2*PI, 0.317188}, {2.763005-2*PI, 0.283019}, \
        {2.799996-2*PI, 0.249716}, {2.838730-2*PI, 0.217321}, {2.879130-2*PI, 0.185875}, \
        {2.921118-2*PI, 0.155417}, {2.964609-2*PI, 0.125982}, {3.009514-2*PI, 0.097608}, \
        {3.055741-2*PI, 0.070327}, {3.103194-2*PI, 0.044172}, {-3.131412, 0.019175}, \
        {3.093238-2*PI, 0.135772}, {-3.134655, 0.104372}, {-3.078635, 0.074316}, \
        {-3.022024, 0.045653}, {-2.964963, 0.018428}, {-2.914903, 0.007316}, \
        {-2.881573, 0.031542}, {-2.844565, 0.050054}, {-2.810218, 0.075297}, \
        {-2.772181, 0.094339}, {-2.733275, 0.112594}, {-2.692963, 0.128759}, \
        {-2.651559, 0.143242}, {-2.609164, 0.156027}, {-2.565879, 0.167102}, \
        {-2.521802, 0.176457}, {-2.477029, 0.184087}, {-2.431657, 0.189988}, \
        {-2.385776, 0.194160}, {-2.339479, 0.196605}, {-2.292850, 0.197328}, \
        {-2.245976, 0.196339}, {-2.198937, 0.193647}, {-2.151810, 0.189266}, \
        {-2.104669, 0.183211}, {-2.057585, 0.175501}, {-2.010625, 0.166155}, \
        {-1.963852, 0.155197}, {-1.917324, 0.142652}, {-1.871098, 0.128544}, \
        {-1.825224, 0.112905}, {-1.779751, 0.095762}, {-1.734723, 0.077150}, \
        {-1.690181, 0.057102}, {-1.646161, 0.035654}, {-1.602698, 0.012842}, \
        {-1.571117, 0.011294}, {-1.554279, 0.036715}, {-1.539324, 0.063379}, \
        {-1.526232, 0.091243}, {-1.514980, 0.120263}, {-1.505541, 0.150394}, \
        {-1.497884, 0.181589}, {-1.491977, 0.213800}, {-1.487784, 0.246979}, \
        {-1.485267, 0.281075}, {-1.484388, 0.316039}, {-1.485105, 0.351817}, \
        {-1.487377, 0.388358}, {-1.491160, 0.425606}, {-1.496410, 0.463509}, \
        {-1.503083, 0.502009}, {-1.511135, 0.541051}, {-1.520521, 0.580577}, \
        {-1.531198, 0.620529}, {-1.543120, 0.660849}, {-1.556246, 0.701477}, \
        {-1.570533, 0.742354}, {-1.585940, 0.783417}, {-1.602428, 0.824607}, \
        {-1.619957, 0.865861}, {-1.638490, 0.907118}, {-1.657993, 0.948314}, \
        {-1.678430, 0.989386}, {-1.699770, 1.030272}, {-1.721981, 1.070909}, \
        {-1.745035, 1.111232}, {-1.768904, 1.151178}, {-1.793563, 1.190684}, \
        {-1.818987, 1.229686}, {-1.845154, 1.268122}, {-1.872042, 1.305928}, \
        {-1.899633, 1.343043}, {-1.927907, 1.379405}, {-1.956847, 1.414953}, \
        {-1.986437, 1.449627}, {-2.016661, 1.483368}, {-2.047505, 1.516117}, \
        {-2.078955, 1.547818}, {-2.110996, 1.578415}, {-2.143615, 1.607854}, \
        {-2.176799, 1.636082}, {-2.210533, 1.663049}, {-2.244803, 1.688705}, \
        {-2.279594, 1.713003}, {-2.314891, 1.735898}, {-2.350677, 1.757349}, \
        {-2.386934, 1.777313}, {-2.423643, 1.795753}, {-2.460784, 1.812634}};

double hindlimb[ROWS][COLS] = {{-0.865582, -0.250746}, {-0.839609, -0.233767}, {-0.813066, -0.219054}, \
        {-0.785986, -0.206649}, {-0.758404, -0.196587}, {-0.730357, -0.188892}, \
        {-0.701890, -0.183580}, {-0.673047, -0.180655}, {-0.643880, -0.180115}, \
        {-0.614441, -0.181945}, {-0.584784, -0.186124}, {-0.554969, -0.192619}, \
        {-0.525056, -0.201389}, {-0.495106, -0.212383}, {-0.465185, -0.225543}, \
        {-0.435356, -0.240801}, {-0.405687, -0.258083}, {-0.376244, -0.277305}, \
        {-0.347094, -0.298378}, {-0.318305, -0.321206}, {-0.289942, -0.345686}, \
        {-0.262071, -0.371710}, {-0.234758, -0.399167}, {-0.208065, -0.427940}, \
        {-0.182054, -0.457908}, {-0.156784, -0.488949}, {-0.132313, -0.520938}, \
        {-0.108694, -0.553749}, {-0.085980, -0.587255}, {-0.064220, -0.621330}, \
        {-0.043457, -0.655846}, {-0.023734, -0.690679}, {-0.005090, -0.725707}, \
        {0.012443, -0.760810}, {0.028834, -0.795872}, {0.044055, -0.830780}, \
        {0.058087, -0.865425}, {0.070909, -0.899706}, {0.082509, -0.933525}, \
        {0.092878, -0.966790}, {0.102008, -0.999418}, {0.109899, -1.031329}, \
        {0.116554, -1.062453}, {0.121978, -1.092727}, {0.126181, -1.122094}, \
        {0.129178, -1.150506}, {0.130984, -1.177922}, {0.131621, -1.204310}, \
        {0.131110, -1.229644}, {0.129478, -1.253908}, {0.126754, -1.277092}, \
        {0.122968, -1.299193}, {0.118152, -1.320217}, {0.112340, -1.340175}, \
        {0.105570, -1.359086}, {0.097877, -1.376974}, {0.089299, -1.393869}, \
        {0.079876, -1.409808}, {0.069646, -1.424830}, {0.058648, -1.438980}, \
        {0.046922, -1.452308}, {0.034506, -1.464865}, {0.021437, -1.476705}, \
        {0.007754, -1.487884}, {-0.006508, -1.498461}, {-0.021314, -1.508494}, \
        {-0.036630, -1.518042}, {-0.052424, -1.527163}, {-0.068666, -1.535916}, \
        {-0.085328, -1.544354}, {-0.102382, -1.552533}, {-0.119805, -1.560503}, \
        {-0.137574, -1.568311}, {-0.155667, -1.576002}, {-0.174067, -1.583616}, \
        {-0.192756, -1.591187}, {-0.211721, -1.598747}, {-0.230948, -1.606319}, \
        {-0.250427, -1.613923}, {-0.270148, -1.621574}, {-0.290104, -1.629277}, \
        {-0.310289, -1.637035}, {-0.330697, -1.644843}, {-0.351325, -1.652690}, \
        {-0.372170, -1.660558}, {-0.393229, -1.668424}, {-0.414502, -1.676257}, \
        {-0.435986, -1.684024}, {-0.457681, -1.691681}, {-0.479584, -1.699184}, \
        {-0.501694, -1.706479}, {-0.524009, -1.713512}, {-0.546525, -1.720221}, \
        {-0.569237, -1.726541}, {-0.592140, -1.732406}, {-0.615226, -1.737745}, \
        {-0.638487, -1.742484}, {-0.661912, -1.746549}, {-0.685487, -1.749864}, \
        {-0.709197, -1.752354}, {-0.733025, -1.753941}, {-0.756952, -1.754551}, \
        {-0.780953, -1.754109}, {-0.805004, -1.752544}, {-0.829077, -1.749785}, \
        {-0.853142, -1.745767}, {-0.877164, -1.740427}, {-0.901107, -1.733707}, \
        {-0.924932, -1.725554}, {-0.948598, -1.715920}, {-0.972060, -1.704763}, \
        {-0.995272, -1.692048}, {-1.018185, -1.677746}, {-1.040747, -1.661837}, \
        {-1.062906, -1.644306}, {-1.084606, -1.625146}, {-1.105791, -1.604362}, \
        {-1.126403, -1.581961}, {-1.146384, -1.557964}, {-1.165675, -1.532397}, \
        {-1.184215, -1.505295}, {-1.201944, -1.476703}, {-1.218802, -1.446673}, \
        {-1.234731, -1.415264}, {-1.249671, -1.382547}, {-1.263565, -1.348595}, \
        {-1.276357, -1.313494}, {-1.287994, -1.277334}, {-1.298424, -1.240211}, \
        {-1.307597, -1.202230}, {-1.315468, -1.163499}, {-1.321992, -1.124134}, \
        {-1.327131, -1.084253}, {-1.330848, -1.043979}, {-1.333111, -1.003438}, \
        {-1.333894, -0.962761}, {-1.333172, -0.922079}, {-1.330928, -0.881524}, \
        {-1.327148, -0.841231}, {-1.321824, -0.801333}, {-1.314953, -0.761963}, \
        {-1.306538, -0.723252}, {-1.296586, -0.685332}, {-1.285110, -0.648329}, \
        {-1.272131, -0.612366}, {-1.257671, -0.577563}, {-1.241761, -0.544036}, \
        {-1.224436, -0.511894}, {-1.205737, -0.481242}, {-1.185710, -0.452177}, \
        {-1.164405, -0.424792}, {-1.141879, -0.399168}, {-1.118192, -0.375384}, \
        {-1.093409, -0.353508}, {-1.067599, -0.333600}, {-1.040834, -0.315712}, \
        {-1.013192, -0.299888}, {-0.984752, -0.286161}, {-0.955597, -0.274557}, \
        {-0.925812, -0.265094}, {-0.895484, -0.257779}, {-0.864703, -0.252611}};

double hindlimb_cont[ROWS][COLS] = {{-0.351325, -1.652690}, {-0.372170, -1.660558}, {-0.393229, -1.668424}, \
        {-0.414502, -1.676257}, {-0.435986, -1.684024}, {-0.457681, -1.691681}, \
        {-0.479584, -1.699184}, {-0.501694, -1.706479}, {-0.524009, -1.713512}, \
        {-0.546525, -1.720221}, {-0.569237, -1.726542}, {-0.592140, -1.732406}, \
        {-0.615226, -1.737745}, {-0.638487, -1.742484}, {-0.661912, -1.746550}, \
        {-0.685487, -1.749864}, {-0.709197, -1.752354}, {-0.733025, -1.753941}, \
        {-0.756952, -1.754551}, {-0.780953, -1.754109}, {-0.805004, -1.752544}, \
        {-0.829077, -1.749785}, {-0.853142, -1.745767}, {-0.877164, -1.740427}, \
        {-0.901107, -1.733707}, {-0.924932, -1.725554}, {-0.948598, -1.715920}, \
        {-0.972060, -1.704763}, {-0.995272, -1.692048}, {-1.018185, -1.677746}, \
        {-1.040747, -1.661837}, {-1.062906, -1.644306}, {-1.084606, -1.625146}, \
        {-1.105791, -1.604362}, {-1.126403, -1.581961}, {-1.146384, -1.557964}, \
        {-1.165675, -1.532397}, {-1.184215, -1.505295}, {-1.201944, -1.476703}, \
        {-1.218802, -1.446673}, {-1.234731, -1.415264}, {-1.249671, -1.382547}, \
        {-1.263565, -1.348595}, {-1.276357, -1.313494}, {-1.287994, -1.277334}, \
        {-1.298424, -1.240211}, {-1.307597, -1.202230}, {-1.315468, -1.163499}, \
        {-1.321992, -1.124134}, {-1.327131, -1.084253}, {-1.330848, -1.043979}, \
        {-1.333111, -1.003438}, {-1.333894, -0.962761}, {-1.333172, -0.922079}, \
        {-1.330928, -0.881524}, {-1.327148, -0.841231}, {-1.321824, -0.801333}, \
        {-1.314953, -0.761963}, {-1.306538, -0.723252}, {-1.296586, -0.685332}, \
        {-1.285110, -0.648329}, {-1.272131, -0.612366}, {-1.257671, -0.577563}, \
        {-1.241761, -0.544036}, {-1.224436, -0.511894}, {-1.205737, -0.481242}, \
        {-1.185710, -0.452177}, {-1.164405, -0.424792}, {-1.141879, -0.399168}, \
        {-1.118192, -0.375384}, {-1.093409, -0.353508}, {-1.067599, -0.333600}, \
        {-1.040834, -0.315712}, {-1.013192, -0.299888}, {-0.984752, -0.286161}, \
        {-0.955597, -0.274557}, {-0.925812, -0.265094}, {-0.895484, -0.257779}, \
        {-0.864703, -0.252611}, {-0.865582, -0.250746}, {-0.839609, -0.233767}, \
        {-0.813066, -0.219054}, {-0.785986, -0.206649}, {-0.758404, -0.196587}, \
        {-0.730357, -0.188892}, {-0.701890, -0.183580}, {-0.673047, -0.180655}, \
        {-0.643880, -0.180115}, {-0.614441, -0.181945}, {-0.584784, -0.186124}, \
        {-0.554969, -0.192619}, {-0.525056, -0.201389}, {-0.495106, -0.212383}, \
        {-0.465185, -0.225543}, {-0.435356, -0.240801}, {-0.405687, -0.258083}, \
        {-0.376244, -0.277305}, {-0.347094, -0.298378}, {-0.318305, -0.321206}, \
        {-0.289942, -0.345686}, {-0.262071, -0.371710}, {-0.234758, -0.399167}, \
        {-0.208065, -0.427940}, {-0.182054, -0.457908}, {-0.156784, -0.488949}, \
        {-0.132313, -0.520938}, {-0.108694, -0.553749}, {-0.085980, -0.587255}, \
        {-0.064220, -0.621330}, {-0.043457, -0.655846}, {-0.023734, -0.690679}, \
        {-0.005090, -0.725707}, {0.012443, -0.760810}, {0.028834, -0.795872}, \
        {0.044055, -0.830780}, {0.058087, -0.865425}, {0.070909, -0.899706}, \
        {0.082509, -0.933525}, {0.092878, -0.966790}, {0.102008, -0.999418}, \
        {0.109899, -1.031329}, {0.116554, -1.062453}, {0.121978, -1.092727}, \
        {0.126181, -1.122094}, {0.129178, -1.150506}, {0.130984, -1.177922}, \
        {0.131621, -1.204311}, {0.131110, -1.229644}, {0.129478, -1.253908}, \
        {0.126754, -1.277092}, {0.122968, -1.299193}, {0.118152, -1.320217}, \
        {0.112340, -1.340175}, {0.105570, -1.359086}, {0.097877, -1.376974}, \
        {0.089299, -1.393869}, {0.079876, -1.409808}, {0.069646, -1.424830}, \
        {0.058648, -1.438980}, {0.046922, -1.452308}, {0.034506, -1.464865}, \
        {0.021437, -1.476705}, {0.007754, -1.487884}, {-0.006508, -1.498461}, \
        {-0.021314, -1.508494}, {-0.036630, -1.518042}, {-0.052424, -1.527163}, \
        {-0.068666, -1.535916}, {-0.085328, -1.544354}, {-0.102382, -1.552533}, \
        {-0.119805, -1.560503}, {-0.137574, -1.568311}, {-0.155667, -1.576002}, \
        {-0.174067, -1.583616}, {-0.192756, -1.591187}, {-0.211721, -1.598747}, \
        {-0.230948, -1.606319}, {-0.250427, -1.613923}, {-0.270148, -1.621574}, \
        {-0.290104, -1.629277}, {-0.310289, -1.637035}, {-0.330697, -1.644843}};


void getGait(double _dst[4][2], double _Ts, double _t, double _start_t)
{
        static int count = 0;
        if (_start_t == 0)
        {
                count = 0;
        }
        else
        {
                double t = (_t - _start_t) - _Ts * int((_t - _start_t) / _Ts);
                double step_t = _Ts / double(ROWS);
                count = int(t / step_t);
                count = (count < 0) ? 0 : ((count > (ROWS - 1)) ? (ROWS - 1) : count);
        }
        _dst[0][0] = forelimb[count][0];
        _dst[0][1] = forelimb[count][1];
        _dst[1][0] = forelimb_cont[count][0];
        _dst[1][1] = forelimb_cont[count][1];
        _dst[2][0] = hindlimb[count][0];
        _dst[2][1] = hindlimb[count][1];
        _dst[3][0] = hindlimb_cont[count][0];
        _dst[3][1] = hindlimb_cont[count][1];
}
