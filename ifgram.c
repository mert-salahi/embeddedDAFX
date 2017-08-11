#include <stdio.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include "ifgram.h"
#include "ifptrack.h"


float wu[N];
#pragma DATA_SECTION(wu, ".EXT_RAM")

float du[N];
#pragma DATA_SECTION(du, ".EXT_RAM")


//Hanning window coefficients
float win[1024] = {0,9.4124e-06,3.7649e-05,8.4709e-05,0.00015059,0.00023529,0.00033881,0.00046114,0.00060227,0.00076221,0.00094094,0.0011385,0.0013548,0.0015899,0.0018437,0.0021163,0.0024076,0.0027177,0.0030465,0.003394,0.0037602,0.0041451,0.0045487,0.0049709,0.0054117,0.0058712,0.0063493,0.006846,0.0073612,0.007895,0.0084473,0.0090181,0.0096074,0.010215,0.010841,0.011486,0.012149,0.01283,0.01353,0.014248,0.014984,0.015739,0.016512,0.017303,0.018112,0.018939,0.019785,0.020648,0.02153,0.022429,0.023347,0.024282,0.025236,0.026207,0.027196,0.028203,0.029228,0.03027,0.03133,0.032408,0.033504,0.034617,0.035747,0.036895,0.03806,0.039243,0.040443,0.04166,0.042895,0.044147,0.045416,0.046702,0.048005,0.049326,0.050663,0.052017,0.053388,0.054776,0.05618,0.057601,0.059039,0.060494,0.061965,0.063453,0.064957,0.066477,0.068014,0.069567,0.071136,0.072721,0.074322,0.07594,0.077573,0.079223,0.080888,0.082569,0.084265,0.085977,0.087705,0.089449,0.091208,0.092982,0.094771,0.096576,0.098396,0.10023,0.10208,0.10395,0.10583,0.10772,0.10963,0.11156,0.11349,0.11545,0.11742,0.1194,0.1214,0.12341,0.12543,0.12747,0.12952,0.13159,0.13367,0.13577,0.13788,0.14,0.14213,0.14428,0.14645,0.14862,0.15081,0.15301,0.15523,0.15746,0.1597,0.16195,0.16422,0.1665,0.16879,0.1711,0.17341,0.17574,0.17808,0.18044,0.1828,0.18518,0.18757,0.18997,0.19238,0.19481,0.19724,0.19969,0.20215,0.20462,0.2071,0.20959,0.2121,0.21461,0.21713,0.21967,0.22221,0.22477,0.22734,0.22991,0.2325,0.2351,0.23771,0.24032,0.24295,0.24558,0.24823,0.25089,0.25355,0.25622,0.25891,0.2616,0.2643,0.26701,0.26973,0.27246,0.27519,0.27794,0.28069,0.28345,0.28622,0.289,0.29179,0.29458,0.29738,0.30019,0.303,0.30583,0.30866,0.3115,0.31434,0.31719,0.32005,0.32292,0.32579,0.32867,0.33156,0.33445,0.33734,0.34025,0.34316,0.34608,0.349,0.35192,0.35486,0.3578,0.36074,0.36369,0.36664,0.3696,0.37257,0.37554,0.37851,0.38149,0.38447,0.38746,0.39045,0.39344,0.39644,0.39945,0.40245,0.40547,0.40848,0.4115,0.41452,0.41754,0.42057,0.4236,0.42663,0.42967,0.43271,0.43575,0.43879,0.44184,0.44489,0.44794,0.45099,0.45405,0.4571,0.46016,0.46322,0.46628,0.46934,0.4724,0.47547,0.47853,0.4816,0.48466,0.48773,0.4908,0.49386,0.49693,0.5,0.50307,0.50614,0.5092,0.51227,0.51534,0.5184,0.52147,0.52453,0.5276,0.53066,0.53372,0.53678,0.53984,0.5429,0.54595,0.54901,0.55206,0.55511,0.55816,0.56121,0.56425,0.56729,0.57033,0.57337,0.5764,0.57943,0.58246,0.58548,0.5885,0.59152,0.59453,0.59755,0.60055,0.60356,0.60656,0.60955,0.61254,0.61553,0.61851,0.62149,0.62446,0.62743,0.6304,0.63336,0.63631,0.63926,0.6422,0.64514,0.64808,0.651,0.65392,0.65684,0.65975,0.66266,0.66555,0.66844,0.67133,0.67421,0.67708,0.67995,0.68281,0.68566,0.6885,0.69134,0.69417,0.697,0.69981,0.70262,0.70542,0.70821,0.711,0.71378,0.71655,0.71931,0.72206,0.72481,0.72754,0.73027,0.73299,0.7357,0.7384,0.74109,0.74378,0.74645,0.74911,0.75177,0.75442,0.75705,0.75968,0.76229,0.7649,0.7675,0.77009,0.77266,0.77523,0.77779,0.78033,0.78287,0.78539,0.7879,0.79041,0.7929,0.79538,0.79785,0.80031,0.80276,0.80519,0.80762,0.81003,0.81243,0.81482,0.8172,0.81956,0.82192,0.82426,0.82659,0.8289,0.83121,0.8335,0.83578,0.83805,0.8403,0.84254,0.84477,0.84699,0.84919,0.85138,0.85355,0.85572,0.85787,0.86,0.86212,0.86423,0.86633,0.86841,0.87048,0.87253,0.87457,0.87659,0.8786,0.8806,0.88258,0.88455,0.88651,0.88844,0.89037,0.89228,0.89417,0.89605,0.89792,0.89977,0.9016,0.90342,0.90523,0.90702,0.90879,0.91055,0.91229,0.91402,0.91573,0.91743,0.91911,0.92078,0.92243,0.92406,0.92568,0.92728,0.92886,0.93043,0.93199,0.93352,0.93504,0.93655,0.93804,0.93951,0.94096,0.9424,0.94382,0.94522,0.94661,0.94798,0.94934,0.95067,0.95199,0.9533,0.95458,0.95585,0.9571,0.95834,0.95956,0.96076,0.96194,0.96311,0.96425,0.96538,0.9665,0.96759,0.96867,0.96973,0.97077,0.9718,0.9728,0.97379,0.97476,0.97572,0.97665,0.97757,0.97847,0.97935,0.98022,0.98106,0.98189,0.9827,0.98349,0.98426,0.98502,0.98575,0.98647,0.98717,0.98785,0.98851,0.98916,0.98978,0.99039,0.99098,0.99155,0.99211,0.99264,0.99315,0.99365,0.99413,0.99459,0.99503,0.99545,0.99585,0.99624,0.99661,0.99695,0.99728,0.99759,0.99788,0.99816,0.99841,0.99865,0.99886,0.99906,0.99924,0.9994,0.99954,0.99966,0.99976,0.99985,0.99992,0.99996,0.99999,1,0.99999,0.99996,0.99992,0.99985,0.99976,0.99966,0.99954,0.9994,0.99924,0.99906,0.99886,0.99865,0.99841,0.99816,0.99788,0.99759,0.99728,0.99695,0.99661,0.99624,0.99585,0.99545,0.99503,0.99459,0.99413,0.99365,0.99315,0.99264,0.99211,0.99155,0.99098,0.99039,0.98978,0.98916,0.98851,0.98785,0.98717,0.98647,0.98575,0.98502,0.98426,0.98349,0.9827,0.98189,0.98106,0.98022,0.97935,0.97847,0.97757,0.97665,0.97572,0.97476,0.97379,0.9728,0.9718,0.97077,0.96973,0.96867,0.96759,0.9665,0.96538,0.96425,0.96311,0.96194,0.96076,0.95956,0.95834,0.9571,0.95585,0.95458,0.9533,0.95199,0.95067,0.94934,0.94798,0.94661,0.94522,0.94382,0.9424,0.94096,0.93951,0.93804,0.93655,0.93504,0.93352,0.93199,0.93043,0.92886,0.92728,0.92568,0.92406,0.92243,0.92078,0.91911,0.91743,0.91573,0.91402,0.91229,0.91055,0.90879,0.90702,0.90523,0.90342,0.9016,0.89977,0.89792,0.89605,0.89417,0.89228,0.89037,0.88844,0.88651,0.88455,0.88258,0.8806,0.8786,0.87659,0.87457,0.87253,0.87048,0.86841,0.86633,0.86423,0.86212,0.86,0.85787,0.85572,0.85355,0.85138,0.84919,0.84699,0.84477,0.84254,0.8403,0.83805,0.83578,0.8335,0.83121,0.8289,0.82659,0.82426,0.82192,0.81956,0.8172,0.81482,0.81243,0.81003,0.80762,0.80519,0.80276,0.80031,0.79785,0.79538,0.7929,0.79041,0.7879,0.78539,0.78287,0.78033,0.77779,0.77523,0.77266,0.77009,0.7675,0.7649,0.76229,0.75968,0.75705,0.75442,0.75177,0.74911,0.74645,0.74378,0.74109,0.7384,0.7357,0.73299,0.73027,0.72754,0.72481,0.72206,0.71931,0.71655,0.71378,0.711,0.70821,0.70542,0.70262,0.69981,0.697,0.69417,0.69134,0.6885,0.68566,0.68281,0.67995,0.67708,0.67421,0.67133,0.66844,0.66555,0.66266,0.65975,0.65684,0.65392,0.651,0.64808,0.64514,0.6422,0.63926,0.63631,0.63336,0.6304,0.62743,0.62446,0.62149,0.61851,0.61553,0.61254,0.60955,0.60656,0.60356,0.60055,0.59755,0.59453,0.59152,0.5885,0.58548,0.58246,0.57943,0.5764,0.57337,0.57033,0.56729,0.56425,0.56121,0.55816,0.55511,0.55206,0.54901,0.54595,0.5429,0.53984,0.53678,0.53372,0.53066,0.5276,0.52453,0.52147,0.5184,0.51534,0.51227,0.5092,0.50614,0.50307,0.5,0.49693,0.49386,0.4908,0.48773,0.48466,0.4816,0.47853,0.47547,0.4724,0.46934,0.46628,0.46322,0.46016,0.4571,0.45405,0.45099,0.44794,0.44489,0.44184,0.43879,0.43575,0.43271,0.42967,0.42663,0.4236,0.42057,0.41754,0.41452,0.4115,0.40848,0.40547,0.40245,0.39945,0.39644,0.39344,0.39045,0.38746,0.38447,0.38149,0.37851,0.37554,0.37257,0.3696,0.36664,0.36369,0.36074,0.3578,0.35486,0.35192,0.349,0.34608,0.34316,0.34025,0.33734,0.33445,0.33156,0.32867,0.32579,0.32292,0.32005,0.31719,0.31434,0.3115,0.30866,0.30583,0.303,0.30019,0.29738,0.29458,0.29179,0.289,0.28622,0.28345,0.28069,0.27794,0.27519,0.27246,0.26973,0.26701,0.2643,0.2616,0.25891,0.25622,0.25355,0.25089,0.24823,0.24558,0.24295,0.24032,0.23771,0.2351,0.2325,0.22991,0.22734,0.22477,0.22221,0.21967,0.21713,0.21461,0.2121,0.20959,0.2071,0.20462,0.20215,0.19969,0.19724,0.19481,0.19238,0.18997,0.18757,0.18518,0.1828,0.18044,0.17808,0.17574,0.17341,0.1711,0.16879,0.1665,0.16422,0.16195,0.1597,0.15746,0.15523,0.15301,0.15081,0.14862,0.14645,0.14428,0.14213,0.14,0.13788,0.13577,0.13367,0.13159,0.12952,0.12747,0.12543,0.12341,0.1214,0.1194,0.11742,0.11545,0.11349,0.11156,0.10963,0.10772,0.10583,0.10395,0.10208,0.10023,0.098396,0.096576,0.094771,0.092982,0.091208,0.089449,0.087705,0.085977,0.084265,0.082569,0.080888,0.079223,0.077573,0.07594,0.074322,0.072721,0.071136,0.069567,0.068014,0.066477,0.064957,0.063453,0.061965,0.060494,0.059039,0.057601,0.05618,0.054776,0.053388,0.052017,0.050663,0.049326,0.048005,0.046702,0.045416,0.044147,0.042895,0.04166,0.040443,0.039243,0.03806,0.036895,0.035747,0.034617,0.033504,0.032408,0.03133,0.03027,0.029228,0.028203,0.027196,0.026207,0.025236,0.024282,0.023347,0.022429,0.02153,0.020648,0.019785,0.018939,0.018112,0.017303,0.016512,0.015739,0.014984,0.014248,0.01353,0.01283,0.012149,0.011486,0.010841,0.010215,0.0096074,0.0090181,0.0084473,0.007895,0.0073612,0.006846,0.0063493,0.0058712,0.0054117,0.0049709,0.0045487,0.0041451,0.0037602,0.003394,0.0030465,0.0027177,0.0024076,0.0021163,0.0018437,0.0015899,0.0013548,0.0011385,0.00094094,0.00076221,0.00060227,0.00046114,0.00033881,0.00023529,0.00015059,8.4709e-05,3.7649e-05,9.4124e-06};
#pragma DATA_SECTION(win, ".EXT_RAM")

//Differentiated hanning window coefficients
float dwin[1024] = {-0,-0.30119,-0.60238,-0.90354,-1.2047,-1.5057,-1.8068,-2.1077,-2.4086,-2.7094,-3.0101,-3.3106,-3.6111,-3.9114,-4.2116,-4.5116,-4.8114,-5.1111,-5.4105,-5.7098,-6.0088,-6.3076,-6.6062,-6.9045,-7.2026,-7.5004,-7.7979,-8.0952,-8.3921,-8.6887,-8.9849,-9.2809,-9.5765,-9.8717,-10.167,-10.461,-10.755,-11.049,-11.342,-11.635,-11.927,-12.219,-12.511,-12.802,-13.092,-13.382,-13.672,-13.961,-14.249,-14.537,-14.825,-15.112,-15.398,-15.684,-15.969,-16.253,-16.537,-16.82,-17.103,-17.385,-17.666,-17.947,-18.227,-18.506,-18.785,-19.063,-19.34,-19.617,-19.892,-20.167,-20.441,-20.715,-20.988,-21.259,-21.531,-21.801,-22.07,-22.339,-22.607,-22.874,-23.14,-23.405,-23.669,-23.933,-24.195,-24.457,-24.717,-24.977,-25.236,-25.494,-25.751,-26.007,-26.262,-26.516,-26.769,-27.021,-27.271,-27.521,-27.77,-28.018,-28.265,-28.511,-28.755,-28.999,-29.241,-29.483,-29.723,-29.962,-30.2,-30.437,-30.673,-30.907,-31.141,-31.373,-31.604,-31.834,-32.063,-32.29,-32.516,-32.741,-32.965,-33.188,-33.409,-33.629,-33.848,-34.065,-34.281,-34.496,-34.71,-34.922,-35.133,-35.343,-35.551,-35.758,-35.964,-36.168,-36.371,-36.573,-36.773,-36.972,-37.169,-37.365,-37.56,-37.753,-37.945,-38.135,-38.324,-38.512,-38.698,-38.882,-39.066,-39.247,-39.427,-39.606,-39.783,-39.959,-40.133,-40.306,-40.477,-40.647,-40.815,-40.981,-41.146,-41.31,-41.472,-41.632,-41.791,-41.948,-42.104,-42.258,-42.41,-42.561,-42.71,-42.858,-43.004,-43.148,-43.291,-43.432,-43.572,-43.71,-43.846,-43.981,-44.114,-44.245,-44.374,-44.502,-44.629,-44.753,-44.876,-44.997,-45.117,-45.235,-45.351,-45.465,-45.578,-45.689,-45.798,-45.906,-46.012,-46.116,-46.218,-46.319,-46.417,-46.514,-46.61,-46.703,-46.795,-46.885,-46.974,-47.06,-47.145,-47.228,-47.309,-47.389,-47.466,-47.542,-47.616,-47.689,-47.759,-47.828,-47.895,-47.96,-48.023,-48.085,-48.144,-48.202,-48.258,-48.312,-48.365,-48.415,-48.464,-48.511,-48.556,-48.599,-48.641,-48.68,-48.718,-48.754,-48.788,-48.821,-48.851,-48.88,-48.906,-48.931,-48.954,-48.976,-48.995,-49.013,-49.028,-49.042,-49.054,-49.064,-49.073,-49.079,-49.084,-49.086,-49.087,-49.086,-49.084,-49.079,-49.073,-49.064,-49.054,-49.042,-49.028,-49.013,-48.995,-48.976,-48.954,-48.931,-48.906,-48.88,-48.851,-48.821,-48.788,-48.754,-48.718,-48.68,-48.641,-48.599,-48.556,-48.511,-48.464,-48.415,-48.365,-48.312,-48.258,-48.202,-48.144,-48.085,-48.023,-47.96,-47.895,-47.828,-47.759,-47.689,-47.616,-47.542,-47.466,-47.389,-47.309,-47.228,-47.145,-47.06,-46.974,-46.885,-46.795,-46.703,-46.61,-46.514,-46.417,-46.319,-46.218,-46.116,-46.012,-45.906,-45.798,-45.689,-45.578,-45.465,-45.351,-45.235,-45.117,-44.997,-44.876,-44.753,-44.629,-44.502,-44.374,-44.245,-44.114,-43.981,-43.846,-43.71,-43.572,-43.432,-43.291,-43.148,-43.004,-42.858,-42.71,-42.561,-42.41,-42.258,-42.104,-41.948,-41.791,-41.632,-41.472,-41.31,-41.146,-40.981,-40.815,-40.647,-40.477,-40.306,-40.133,-39.959,-39.783,-39.606,-39.427,-39.247,-39.066,-38.882,-38.698,-38.512,-38.324,-38.135,-37.945,-37.753,-37.56,-37.365,-37.169,-36.972,-36.773,-36.573,-36.371,-36.168,-35.964,-35.758,-35.551,-35.343,-35.133,-34.922,-34.71,-34.496,-34.281,-34.065,-33.848,-33.629,-33.409,-33.188,-32.965,-32.741,-32.516,-32.29,-32.063,-31.834,-31.604,-31.373,-31.141,-30.907,-30.673,-30.437,-30.2,-29.962,-29.723,-29.483,-29.241,-28.999,-28.755,-28.511,-28.265,-28.018,-27.77,-27.521,-27.271,-27.021,-26.769,-26.516,-26.262,-26.007,-25.751,-25.494,-25.236,-24.977,-24.717,-24.457,-24.195,-23.933,-23.669,-23.405,-23.14,-22.874,-22.607,-22.339,-22.07,-21.801,-21.531,-21.259,-20.988,-20.715,-20.441,-20.167,-19.892,-19.617,-19.34,-19.063,-18.785,-18.506,-18.227,-17.947,-17.666,-17.385,-17.103,-16.82,-16.537,-16.253,-15.969,-15.684,-15.398,-15.112,-14.825,-14.537,-14.249,-13.961,-13.672,-13.382,-13.092,-12.802,-12.511,-12.219,-11.927,-11.635,-11.342,-11.049,-10.755,-10.461,-10.167,-9.8717,-9.5765,-9.2809,-8.9849,-8.6887,-8.3921,-8.0952,-7.7979,-7.5004,-7.2026,-6.9045,-6.6062,-6.3076,-6.0088,-5.7098,-5.4105,-5.1111,-4.8114,-4.5116,-4.2116,-3.9114,-3.6111,-3.3106,-3.0101,-2.7094,-2.4086,-2.1077,-1.8068,-1.5057,-1.2047,-0.90354,-0.60238,-0.30119,-6.0115e-15,0.30119,0.60238,0.90354,1.2047,1.5057,1.8068,2.1077,2.4086,2.7094,3.0101,3.3106,3.6111,3.9114,4.2116,4.5116,4.8114,5.1111,5.4105,5.7098,6.0088,6.3076,6.6062,6.9045,7.2026,7.5004,7.7979,8.0952,8.3921,8.6887,8.9849,9.2809,9.5765,9.8717,10.167,10.461,10.755,11.049,11.342,11.635,11.927,12.219,12.511,12.802,13.092,13.382,13.672,13.961,14.249,14.537,14.825,15.112,15.398,15.684,15.969,16.253,16.537,16.82,17.103,17.385,17.666,17.947,18.227,18.506,18.785,19.063,19.34,19.617,19.892,20.167,20.441,20.715,20.988,21.259,21.531,21.801,22.07,22.339,22.607,22.874,23.14,23.405,23.669,23.933,24.195,24.457,24.717,24.977,25.236,25.494,25.751,26.007,26.262,26.516,26.769,27.021,27.271,27.521,27.77,28.018,28.265,28.511,28.755,28.999,29.241,29.483,29.723,29.962,30.2,30.437,30.673,30.907,31.141,31.373,31.604,31.834,32.063,32.29,32.516,32.741,32.965,33.188,33.409,33.629,33.848,34.065,34.281,34.496,34.71,34.922,35.133,35.343,35.551,35.758,35.964,36.168,36.371,36.573,36.773,36.972,37.169,37.365,37.56,37.753,37.945,38.135,38.324,38.512,38.698,38.882,39.066,39.247,39.427,39.606,39.783,39.959,40.133,40.306,40.477,40.647,40.815,40.981,41.146,41.31,41.472,41.632,41.791,41.948,42.104,42.258,42.41,42.561,42.71,42.858,43.004,43.148,43.291,43.432,43.572,43.71,43.846,43.981,44.114,44.245,44.374,44.502,44.629,44.753,44.876,44.997,45.117,45.235,45.351,45.465,45.578,45.689,45.798,45.906,46.012,46.116,46.218,46.319,46.417,46.514,46.61,46.703,46.795,46.885,46.974,47.06,47.145,47.228,47.309,47.389,47.466,47.542,47.616,47.689,47.759,47.828,47.895,47.96,48.023,48.085,48.144,48.202,48.258,48.312,48.365,48.415,48.464,48.511,48.556,48.599,48.641,48.68,48.718,48.754,48.788,48.821,48.851,48.88,48.906,48.931,48.954,48.976,48.995,49.013,49.028,49.042,49.054,49.064,49.073,49.079,49.084,49.086,49.087,49.086,49.084,49.079,49.073,49.064,49.054,49.042,49.028,49.013,48.995,48.976,48.954,48.931,48.906,48.88,48.851,48.821,48.788,48.754,48.718,48.68,48.641,48.599,48.556,48.511,48.464,48.415,48.365,48.312,48.258,48.202,48.144,48.085,48.023,47.96,47.895,47.828,47.759,47.689,47.616,47.542,47.466,47.389,47.309,47.228,47.145,47.06,46.974,46.885,46.795,46.703,46.61,46.514,46.417,46.319,46.218,46.116,46.012,45.906,45.798,45.689,45.578,45.465,45.351,45.235,45.117,44.997,44.876,44.753,44.629,44.502,44.374,44.245,44.114,43.981,43.846,43.71,43.572,43.432,43.291,43.148,43.004,42.858,42.71,42.561,42.41,42.258,42.104,41.948,41.791,41.632,41.472,41.31,41.146,40.981,40.815,40.647,40.477,40.306,40.133,39.959,39.783,39.606,39.427,39.247,39.066,38.882,38.698,38.512,38.324,38.135,37.945,37.753,37.56,37.365,37.169,36.972,36.773,36.573,36.371,36.168,35.964,35.758,35.551,35.343,35.133,34.922,34.71,34.496,34.281,34.065,33.848,33.629,33.409,33.188,32.965,32.741,32.516,32.29,32.063,31.834,31.604,31.373,31.141,30.907,30.673,30.437,30.2,29.962,29.723,29.483,29.241,28.999,28.755,28.511,28.265,28.018,27.77,27.521,27.271,27.021,26.769,26.516,26.262,26.007,25.751,25.494,25.236,24.977,24.717,24.457,24.195,23.933,23.669,23.405,23.14,22.874,22.607,22.339,22.07,21.801,21.531,21.259,20.988,20.715,20.441,20.167,19.892,19.617,19.34,19.063,18.785,18.506,18.227,17.947,17.666,17.385,17.103,16.82,16.537,16.253,15.969,15.684,15.398,15.112,14.825,14.537,14.249,13.961,13.672,13.382,13.092,12.802,12.511,12.219,11.927,11.635,11.342,11.049,10.755,10.461,10.167,9.8717,9.5765,9.2809,8.9849,8.6887,8.3921,8.0952,7.7979,7.5004,7.2026,6.9045,6.6062,6.3076,6.0088,5.7098,5.4105,5.1111,4.8114,4.5116,4.2116,3.9114,3.6111,3.3106,3.0101,2.7094,2.4086,2.1077,1.8068,1.5057,1.2047,0.90354,0.60238,0.30119};
#pragma DATA_SECTION(dwin, ".EXT_RAM")

//Angular frequency vector
float ww[N] = {0,49.087,98.175,147.26,196.35,245.44,294.52,343.61,392.7,441.79,490.87,539.96,589.05,638.14,687.22,736.31,785.4,834.49,883.57,932.66,981.75,1030.8,1079.9,1129,1178.1,1227.2,1276.3,1325.4,1374.4,1423.5,1472.6,1521.7,1570.8,1619.9,1669,1718.1,1767.1,1816.2,1865.3,1914.4,1963.5,2012.6,2061.7,2110.8,2159.8,2208.9,2258,2307.1,2356.2,2405.3,2454.4,2503.5,2552.5,2601.6,2650.7,2699.8,2748.9,2798,2847.1,2896.2,2945.2,2994.3,3043.4,3092.5,3141.6,3190.7,3239.8,3288.9,3337.9,3387,3436.1,3485.2,3534.3,3583.4,3632.5,3681.6,3730.6,3779.7,3828.8,3877.9,3927,3976.1,4025.2,4074.3,4123.3,4172.4,4221.5,4270.6,4319.7,4368.8,4417.9,4467,4516,4565.1,4614.2,4663.3,4712.4,4761.5,4810.6,4859.7,4908.7,4957.8,5006.9,5056,5105.1,5154.2,5203.3,5252.4,5301.4,5350.5,5399.6,5448.7,5497.8,5546.9,5596,5645,5694.1,5743.2,5792.3,5841.4,5890.5,5939.6,5988.7,6037.7,6086.8,6135.9,6185,6234.1,6283.2,6332.3,6381.4,6430.4,6479.5,6528.6,6577.7,6626.8,6675.9,6725,6774.1,6823.1,6872.2,6921.3,6970.4,7019.5,7068.6,7117.7,7166.8,7215.8,7264.9,7314,7363.1,7412.2,7461.3,7510.4,7559.5,7608.5,7657.6,7706.7,7755.8,7804.9,7854,7903.1,7952.2,8001.2,8050.3,8099.4,8148.5,8197.6,8246.7,8295.8,8344.9,8393.9,8443,8492.1,8541.2,8590.3,8639.4,8688.5,8737.6,8786.6,8835.7,8884.8,8933.9,8983,9032.1,9081.2,9130.3,9179.3,9228.4,9277.5,9326.6,9375.7,9424.8,9473.9,9523,9572,9621.1,9670.2,9719.3,9768.4,9817.5,9866.6,9915.7,9964.7,10014,10063,10112,10161,10210,10259,10308,10357,10407,10456,10505,10554,10603,10652,10701,10750,10799,10848,10897,10946,10996,11045,11094,11143,11192,11241,11290,11339,11388,11437,11486,11536,11585,11634,11683,11732,11781,11830,11879,11928,11977,12026,12075,12125,12174,12223,12272,12321,12370,12419,12468,12517,12566,12615,12665,12714,12763,12812,12861,12910,12959,13008,13057,13106,13155,13205,13254,13303,13352,13401,13450,13499,13548,13597,13646,13695,13744,13794,13843,13892,13941,13990,14039,14088,14137,14186,14235,14284,14334,14383,14432,14481,14530,14579,14628,14677,14726,14775,14824,14873,14923,14972,15021,15070,15119,15168,15217,15266,15315,15364,15413,15463,15512,15561,15610,15659,15708,15757,15806,15855,15904,15953,16002,16052,16101,16150,16199,16248,16297,16346,16395,16444,16493,16542,16592,16641,16690,16739,16788,16837,16886,16935,16984,17033,17082,17131,17181,17230,17279,17328,17377,17426,17475,17524,17573,17622,17671,17721,17770,17819,17868,17917,17966,18015,18064,18113,18162,18211,18261,18310,18359,18408,18457,18506,18555,18604,18653,18702,18751,18800,18850,18899,18948,18997,19046,19095,19144,19193,19242,19291,19340,19390,19439,19488,19537,19586,19635,19684,19733,19782,19831,19880,19929,19979,20028,20077,20126,20175,20224,20273,20322,20371,20420,20469,20519,20568,20617,20666,20715,20764,20813,20862,20911,20960,21009,21058,21108,21157,21206,21255,21304,21353,21402,21451,21500,21549,21598,21648,21697,21746,21795,21844,21893,21942,21991,22040,22089,22138,22187,22237,22286,22335,22384,22433,22482,22531,22580,22629,22678,22727,22777,22826,22875,22924,22973,23022,23071,23120,23169,23218,23267,23317,23366,23415,23464,23513,23562,23611,23660,23709,23758,23807,23856,23906,23955,24004,24053,24102,24151,24200,24249,24298,24347,24396,24446,24495,24544,24593,24642,24691,24740,24789,24838,24887,24936,24985,25035,25084,25133,25182,25231,25280,25329,25378,25427,25476,25525,25575,25624,25673,25722,25771,25820,25869,25918,25967,26016,26065,26114,26164,26213,26262,26311,26360,26409,26458,26507,26556,26605,26654,26704,26753,26802,26851,26900,26949,26998,27047,27096,27145,27194,27243,27293,27342,27391,27440,27489,27538,27587,27636,27685,27734,27783,27833,27882,27931,27980,28029,28078,28127,28176,28225,28274,28323,28373,28422,28471,28520,28569,28618,28667,28716,28765,28814,28863,28912,28962,29011,29060,29109,29158,29207,29256,29305,29354,29403,29452,29502,29551,29600,29649,29698,29747,29796,29845,29894,29943,29992,30041,30091,30140,30189,30238,30287,30336,30385,30434,30483,30532,30581,30631,30680,30729,30778,30827,30876,30925,30974,31023,31072,31121,31170,31220,31269,31318,31367,31416,31465,31514,31563,31612,31661,31710,31760,31809,31858,31907,31956,32005,32054,32103,32152,32201,32250,32299,32349,32398,32447,32496,32545,32594,32643,32692,32741,32790,32839,32889,32938,32987,33036,33085,33134,33183,33232,33281,33330,33379,33429,33478,33527,33576,33625,33674,33723,33772,33821,33870,33919,33968,34018,34067,34116,34165,34214,34263,34312,34361,34410,34459,34508,34558,34607,34656,34705,34754,34803,34852,34901,34950,34999,35048,35097,35147,35196,35245,35294,35343,35392,35441,35490,35539,35588,35637,35687,35736,35785,35834,35883,35932,35981,36030,36079,36128,36177,36226,36276,36325,36374,36423,36472,36521,36570,36619,36668,36717,36766,36816,36865,36914,36963,37012,37061,37110,37159,37208,37257,37306,37356,37405,37454,37503,37552,37601,37650,37699,37748,37797,37846,37895,37945,37994,38043,38092,38141,38190,38239,38288,38337,38386,38435,38485,38534,38583,38632,38681,38730,38779,38828,38877,38926,38975,39024,39074,39123,39172,39221,39270,39319,39368,39417,39466,39515,39564,39614,39663,39712,39761,39810,39859,39908,39957,40006,40055,40104,40153,40203,40252,40301,40350,40399,40448,40497,40546,40595,40644,40693,40743,40792,40841,40890,40939,40988,41037,41086,41135,41184,41233,41282,41332,41381,41430,41479,41528,41577,41626,41675,41724,41773,41822,41872,41921,41970,42019,42068,42117,42166,42215,42264,42313,42362,42412,42461,42510,42559,42608,42657,42706,42755,42804,42853,42902,42951,43001,43050,43099,43148,43197,43246,43295,43344,43393,43442,43491,43541,43590,43639,43688,43737,43786,43835,43884,43933,43982,44031,44080,44130,44179,44228,44277,44326,44375,44424,44473,44522,44571,44620,44670,44719,44768,44817,44866,44915,44964,45013,45062,45111,45160,45209,45259,45308,45357,45406,45455,45504,45553,45602,45651,45700,45749,45799,45848,45897,45946,45995,46044,46093,46142,46191,46240,46289,46338,46388,46437,46486,46535,46584,46633,46682,46731,46780,46829,46878,46928,46977,47026,47075,47124,47173,47222,47271,47320,47369,47418,47468,47517,47566,47615,47664,47713,47762,47811,47860,47909,47958,48007,48057,48106,48155,48204,48253,48302,48351,48400,48449,48498,48547,48597,48646,48695,48744,48793,48842,48891,48940,48989,49038,49087,49136,49186,49235,49284,49333,49382,49431,49480,49529,49578,49627,49676,49726,49775,49824,49873,49922,49971,50020,50069,50118,50167,50216,50265,50315,50364,50413,50462,50511,50560,50609,50658,50707,50756,50805,50855,50904,50953,51002,51051,51100,51149,51198,51247,51296,51345,51394,51444,51493,51542,51591,51640,51689,51738,51787,51836,51885,51934,51984,52033,52082,52131,52180,52229,52278,52327,52376,52425,52474,52524,52573,52622,52671,52720,52769,52818,52867,52916,52965,53014,53063,53113,53162,53211,53260,53309,53358,53407,53456,53505,53554,53603,53653,53702,53751,53800,53849,53898,53947,53996,54045,54094,54143,54192,54242,54291,54340,54389,54438,54487,54536,54585,54634,54683,54732,54782,54831,54880,54929,54978,55027,55076,55125,55174,55223,55272,55321,55371,55420,55469,55518,55567,55616,55665,55714,55763,55812,55861,55911,55960,56009,56058,56107,56156,56205,56254,56303,56352,56401,56450,56500,56549,56598,56647,56696,56745,56794,56843,56892,56941,56990,57040,57089,57138,57187,57236,57285,57334,57383,57432,57481,57530,57580,57629,57678,57727,57776,57825,57874,57923,57972,58021,58070,58119,58169,58218,58267,58316,58365,58414,58463,58512,58561,58610,58659,58709,58758,58807,58856,58905,58954,59003,59052,59101,59150,59199,59248,59298,59347,59396,59445,59494,59543,59592,59641,59690,59739,59788,59838,59887,59936,59985,60034,60083,60132,60181,60230,60279,60328,60377,60427,60476,60525,60574,60623,60672,60721,60770,60819,60868,60917,60967,61016,61065,61114,61163,61212,61261,61310,61359,61408,61457,61506,61556,61605,61654,61703,61752,61801,61850,61899,61948,61997,62046,62096,62145,62194,62243,62292,62341,62390,62439,62488,62537,62586,62636,62685,62734,62783,62832,62881,62930,62979,63028,63077,63126,63175,63225,63274,63323,63372,63421,63470,63519,63568,63617,63666,63715,63765,63814,63863,63912,63961,64010,64059,64108,64157,64206,64255,64304,64354,64403,64452,64501,64550,64599,64648,64697,64746,64795,64844,64894,64943,64992,65041,65090,65139,65188,65237,65286,65335,65384,65433,65483,65532,65581,65630,65679,65728,65777,65826,65875,65924,65973,66023,66072,66121,66170,66219,66268,66317,66366,66415,66464,66513,66562,66612,66661,66710,66759,66808,66857,66906,66955,67004,67053,67102,67152,67201,67250,67299,67348,67397,67446,67495,67544,67593,67642,67692,67741,67790,67839,67888,67937,67986,68035,68084,68133,68182,68231,68281,68330,68379,68428,68477,68526,68575,68624,68673,68722,68771,68821,68870,68919,68968,69017,69066,69115,69164,69213,69262,69311,69360,69410,69459,69508,69557,69606,69655,69704,69753,69802,69851,69900,69950,69999,70048,70097,70146,70195,70244,70293,70342,70391,70440,70489,70539,70588,70637,70686,70735,70784,70833,70882,70931,70980,71029,71079,71128,71177,71226,71275,71324,71373,71422,71471,71520,71569,71618,71668,71717,71766,71815,71864,71913,71962,72011,72060,72109,72158,72208,72257,72306,72355,72404,72453,72502,72551,72600,72649,72698,72748,72797,72846,72895,72944,72993,73042,73091,73140,73189,73238,73287,73337,73386,73435,73484,73533,73582,73631,73680,73729,73778,73827,73877,73926,73975,74024,74073,74122,74171,74220,74269,74318,74367,74416,74466,74515,74564,74613,74662,74711,74760,74809,74858,74907,74956,75006,75055,75104,75153,75202,75251,75300,75349,75398,75447,75496,75545,75595,75644,75693,75742,75791,75840,75889,75938,75987,76036,76085,76135,76184,76233,76282,76331,76380,76429,76478,76527,76576,76625,76674,76724,76773,76822,76871,76920,76969,77018,77067,77116,77165,77214,77264,77313,77362,77411,77460,77509,77558,77607,77656,77705,77754,77804,77853,77902,77951,78000,78049,78098,78147,78196,78245,78294,78343,78393,78442,78491,78540,78589,78638,78687,78736,78785,78834,78883,78933,78982,79031,79080,79129,79178,79227,79276,79325,79374,79423,79472,79522,79571,79620,79669,79718,79767,79816,79865,79914,79963,80012,80062,80111,80160,80209,80258,80307,80356,80405,80454,80503,80552,80601,80651,80700,80749,80798,80847,80896,80945,80994,81043,81092,81141,81191,81240,81289,81338,81387,81436,81485,81534,81583,81632,81681,81730,81780,81829,81878,81927,81976,82025,82074,82123,82172,82221,82270,82320,82369,82418,82467,82516,82565,82614,82663,82712,82761,82810,82860,82909,82958,83007,83056,83105,83154,83203,83252,83301,83350,83399,83449,83498,83547,83596,83645,83694,83743,83792,83841,83890,83939,83989,84038,84087,84136,84185,84234,84283,84332,84381,84430,84479,84528,84578,84627,84676,84725,84774,84823,84872,84921,84970,85019,85068,85118,85167,85216,85265,85314,85363,85412,85461,85510,85559,85608,85657,85707,85756,85805,85854,85903,85952,86001,86050,86099,86148,86197,86247,86296,86345,86394,86443,86492,86541,86590,86639,86688,86737,86786,86836,86885,86934,86983,87032,87081,87130,87179,87228,87277,87326,87376,87425,87474,87523,87572,87621,87670,87719,87768,87817,87866,87916,87965,88014,88063,88112,88161,88210,88259,88308,88357,88406,88455,88505,88554,88603,88652,88701,88750,88799,88848,88897,88946,88995,89045,89094,89143,89192,89241,89290,89339,89388,89437,89486,89535,89584,89634,89683,89732,89781,89830,89879,89928,89977,90026,90075,90124,90174,90223,90272,90321,90370,90419,90468,90517,90566,90615,90664,90713,90763,90812,90861,90910,90959,91008,91057,91106,91155,91204,91253,91303,91352,91401,91450,91499,91548,91597,91646,91695,91744,91793,91842,91892,91941,91990,92039,92088,92137,92186,92235,92284,92333,92382,92432,92481,92530,92579,92628,92677,92726,92775,92824,92873,92922,92972,93021,93070,93119,93168,93217,93266,93315,93364,93413,93462,93511,93561,93610,93659,93708,93757,93806,93855,93904,93953,94002,94051,94101,94150,94199,94248,94297,94346,94395,94444,94493,94542,94591,94640,94690,94739,94788,94837,94886,94935,94984,95033,95082,95131,95180,95230,95279,95328,95377,95426,95475,95524,95573,95622,95671,95720,95769,95819,95868,95917,95966,96015,96064,96113,96162,96211,96260,96309,96359,96408,96457,96506,96555,96604,96653,96702,96751,96800,96849,96898,96948,96997,97046,97095,97144,97193,97242,97291,97340,97389,97438,97488,97537,97586,97635,97684,97733,97782,97831,97880,97929,97978,98028,98077,98126,98175,98224,98273,98322,98371,98420,98469,98518,98567,98617,98666,98715,98764,98813,98862,98911,98960,99009,99058,99107,99157,99206,99255,99304,99353,99402,99451,99500,99549,99598,99647,99696,99746,99795,99844,99893,99942,99991,1.0004e+05,1.0009e+05,1.0014e+05,1.0019e+05,1.0024e+05,1.0029e+05,1.0033e+05,1.0038e+05,1.0043e+05,1.0048e+05};
#pragma DATA_SECTION(ww, ".EXT_RAM")

float fft_du[4097];
#pragma DATA_SECTION(fft_du, ".IRAM")

float fft_wu[4097];
#pragma DATA_SECTION(fft_wu , ".IRAM")


//Sum all the elements in an array
float sumArray(float *a, int numElements){
	int i;
    float sum=0;
    for (i=0; i<numElements; i++)
    {
        sum = sum + *(a+i);
    }
    return sum;

}

//Swap function
void swap(float *v1, float *v2)
{
    float tmp = *v1;
    *v1 = *v2;
    *v2 = tmp;
}


//Function that shifts the zero-frequency component to the center of the matrix
void fftshift(float *data, int count)
{
    int k = 0;
    int c = (int) floor((float)count/2);
    // For odd and for even numbers of element use different algorithm
    if (count % 2 == 0)
    {
        for (k = 0; k < c; k++)
            swap(&data[k], &data[k+c]);
    }
    else
    {
        float tmp = data[0];
        for (k = 0; k < c; k++)
        {
            data[k] = data[c + k + 1];
            data[c + k + 1] = data[k + 1];
        }
        data[c] = tmp;
    }
}


//FFT function
void four1(float data[], int nn, int isign)
{
    int n, mmax, m, j, istep, i;
    float wtemp, wr, wpr, wpi, wi, theta;
    float tempr, tempi;

    n = nn << 1;
    j = 1;
    for (i = 1; i < n; i += 2) {
        if (j > i) {
            tempr = data[j];     data[j] = data[i];     data[i] = tempr;
            tempr = data[j+1]; data[j+1] = data[i+1]; data[i+1] = tempr;
        }
        m = n >> 1;
        while (m >= 2 && j > m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }
    mmax = 2;
    while (n > mmax) {
        istep = 2*mmax;

        theta = 6.283185/(isign*mmax);


        wtemp = sin(0.5*theta);
        wpr = -2.0*wtemp*wtemp;
        wpi = sin(theta);
        wr = 1.0;
        wi = 0.0;
        for (m = 1; m < mmax; m += 2) {
            for (i = m; i <= n; i += istep) {
                j =i + mmax;
                tempr = wr*data[j]   - wi*data[j+1];
                tempi = wr*data[j+1] + wi*data[j];
                data[j]   = data[i]   - tempr;
                data[j+1] = data[i+1] - tempi;
                data[i] += tempr;
                data[i+1] += tempi;
            }
            wr = (wtemp = wr)*wpr - wi*wpi + wr;
            wi = wi*wpr + wtemp*wpi + wi;
        }
        mmax = istep;
    }
}


//DSK board optimized FFT function
void cfftr2_dit(float* x, float* w, short n)
		{
		   short n2, ie, ia, i, j, k, m;
		   float rtemp, itemp, c, s;

		   n2 = n;
		   ie = 1;

		   for(k=n; k > 1; k >>= 1)
		   {
		      n2 >>= 1;
		      ia = 0;
		      for(j=0; j < ie; j++)
		      {
		         c = w[2*j];
		         s = w[2*j+1];
		         for(i=0; i < n2; i++)
		         {
		            m = ia + n2;
		            rtemp     = c * x[2*m]   + s * x[2*m+1];
		            itemp     = c * x[2*m+1] - s * x[2*m];
		            x[2*m]    = x[2*ia]   - rtemp;
		            x[2*m+1]  = x[2*ia+1] - itemp;
		            x[2*ia]   = x[2*ia]   + rtemp;
		            x[2*ia+1] = x[2*ia+1] + itemp;
		            ia++;
		         }
		         ia += n2;
		      }
		      ie <<= 1;
		   }
		}


//ifgram - computes intantaneous frequency and DFT matrices
void ifgram(float *X, float IF[][nhops], float DFT[][nhops]) {
    int i;

    //Look at the portion of signal to prepare for multiplying by W-point window
    //Zeropadding first
    for (i = 0; i<N; i++){
        if (i < nmw1 || i > N-nmw2-1) {
            wu[i]=0;
            du[i]=0;
        }
        else {
            wu[i] = X[i-nmw1]*win[i-nmw1]; //multiply signal by window
            du[i] = X[i-nmw1]*dwin[i-nmw1]; //multiply signal by differentiated window
        }
    }


    int idx1, idx2;
    fft_wu[0] = 0;
    fft_du[0] = 0;

    //Formatting the matrix to be passed into FFT function - Copy the outputs from the windowing into this new matrix. Real components in the odd elements of the matrix and imaginary components in the even elements of the matrix

    for(i=0; i<N/4; i++) {
    	idx1 = (i<<3);
    	idx2 = (i<<2);

        //Unrolling
        fft_wu[idx1+1] = wu[idx2];
        fft_wu[idx1+2] = 0;
        fft_wu[idx1+3] = wu[idx2+1];
        fft_wu[idx1+4] = 0;
        fft_wu[idx1+5] = wu[idx2+2];
        fft_wu[idx1+6] = 0;
        fft_wu[idx1+7] = wu[idx2+3];
        fft_wu[idx1+8] = 0;

        fft_du[idx1+1] = du[idx2];
        fft_du[idx1+2] = 0;
        fft_du[idx1+3] = du[idx2+1];
        fft_du[idx1+4] = 0;
        fft_du[idx1+5] = du[idx2+2];
        fft_du[idx1+6] = 0;
        fft_du[idx1+7] = du[idx2+3];
        fft_du[idx1+8] = 0;
    }

//    /* calculate FFT */

    four1(fft_wu, N, 1);
    four1(fft_du, N, 1);

//    //Scale down to factor out length and window effects
//    //Find absolute values of the complex numbers of fft output and store in DFT array
    for(i=0; i<(N/2+1); i++)
    {
        //Compute absolute value
        float result = sqrtf((fft_wu[2*i+1]*fft_wu[2*i+1])+(fft_wu[2*i+2]*fft_wu[2*i+2]));
        DFT[i][0] = result*norm;
    }

    for(i=0; i<(N/2+1); i++) {
        //Computing the instantaneous frequency from the time derivative of the complex phase spectrum
        float num = (fft_wu[(2*i)+1]*(-fft_du[2*i+2]+ww[i]*fft_wu[2*i+1]))+(fft_wu[(2*i)+2]*(fft_du[2*i+1]+ww[i]*fft_wu[2*i+2]));
        float mag = (fft_wu[(2*i)+1]*fft_wu[(2*i)+1])+(fft_wu[(2*i)+2]*fft_wu[(2*i)+2]);
        float abs = sqrtf(mag);
        bool flag = 0;
        if (abs < 0.00001) {
            flag = 1;
        }

        IF[i][0] = (1/(TWOPI))*num/(mag+flag);
    }
}
