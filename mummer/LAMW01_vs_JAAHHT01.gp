set terminal postscript color solid "Courier" 8
set output "LAMW01_vs_JAAHHT01.ps"
set xtics rotate ( \
 "*gi|816976904|gb|LAMW01000146.1|" 1, \
 "gi|816977202|gb|LAMW01000047.1|" 10363, \
 "*gi|816977288|gb|LAMW01000018.1|" 49579, \
 "gi|816977196|gb|LAMW01000049.1|" 81325, \
 "*gi|816976812|gb|LAMW01000177.1|" 96994, \
 "*gi|816976875|gb|LAMW01000156.1|" 115843, \
 "*gi|816977164|gb|LAMW01000060.1|" 139043, \
 "*gi|816977235|gb|LAMW01000036.1|" 172559, \
 "gi|816977176|gb|LAMW01000056.1|" 184681, \
 "gi|816977105|gb|LAMW01000079.1|" 196823, \
 "*gi|816977220|gb|LAMW01000041.1|" 207566, \
 "gi|816976958|gb|LAMW01000128.1|" 218363, \
 "*gi|816976910|gb|LAMW01000144.1|" 237394, \
 "*gi|816977217|gb|LAMW01000042.1|" 254901, \
 "gi|816977333|gb|LAMW01000003.1|" 280922, \
 "gi|816976997|gb|LAMW01000115.1|" 294340, \
 "*gi|816976746|gb|LAMW01000199.1|" 320554, \
 "gi|816977179|gb|LAMW01000055.1|" 333256, \
 "*gi|816977003|gb|LAMW01000113.1|" 344196, \
 "gi|816977208|gb|LAMW01000045.1|" 376861, \
 "gi|816977276|gb|LAMW01000022.1|" 403484, \
 "*gi|816977158|gb|LAMW01000062.1|" 440283, \
 "*gi|816976976|gb|LAMW01000122.1|" 458705, \
 "*gi|816976767|gb|LAMW01000192.1|" 472315, \
 "gi|816976821|gb|LAMW01000174.1|" 486703, \
 "gi|816977173|gb|LAMW01000057.1|" 501564, \
 "*gi|816976913|gb|LAMW01000143.1|" 519727, \
 "gi|816977099|gb|LAMW01000081.1|" 537474, \
 "*gi|816976928|gb|LAMW01000138.1|" 549698, \
 "gi|816977030|gb|LAMW01000104.1|" 568228, \
 "gi|816977121|gb|LAMW01000073.1|" 580812, \
 "gi|816976881|gb|LAMW01000154.1|" 592086, \
 "*gi|816977167|gb|LAMW01000059.1|" 605160, \
 "gi|816977312|gb|LAMW01000010.1|" 617591, \
 "*gi|816977148|gb|LAMW01000065.1|" 632861, \
 "*gi|816976988|gb|LAMW01000118.1|" 652756, \
 "*gi|816976961|gb|LAMW01000127.1|" 671212, \
 "*gi|816977238|gb|LAMW01000035.1|" 684145, \
 "*gi|816977054|gb|LAMW01000096.1|" 695021, \
 "*gi|816977256|gb|LAMW01000029.1|" 725598, \
 "*gi|816976863|gb|LAMW01000160.1|" 737345, \
 "*gi|816976842|gb|LAMW01000167.1|" 748495, \
 "*gi|816977241|gb|LAMW01000034.1|" 758834, \
 "*gi|816977232|gb|LAMW01000037.1|" 769162, \
 "gi|816976878|gb|LAMW01000155.1|" 797141, \
 "gi|816976916|gb|LAMW01000142.1|" 820257, \
 "gi|816977042|gb|LAMW01000100.1|" 852709, \
 "*gi|816976785|gb|LAMW01000186.1|" 867417, \
 "*gi|816976898|gb|LAMW01000148.1|" 879470, \
 "gi|816976788|gb|LAMW01000185.1|" 896607, \
 "gi|816976752|gb|LAMW01000197.1|" 910407, \
 "*gi|816976979|gb|LAMW01000121.1|" 929900, \
 "*gi|816976848|gb|LAMW01000165.1|" 940100, \
 "*gi|816977119|gb|LAMW01000074.1|" 975791, \
 "*gi|816976779|gb|LAMW01000188.1|" 988589, \
 "gi|816976869|gb|LAMW01000158.1|" 1005686, \
 "*gi|816977253|gb|LAMW01000030.1|" 1017719, \
 "gi|816977132|gb|LAMW01000070.1|" 1030577, \
 "gi|816977111|gb|LAMW01000077.1|" 1044365, \
 "*gi|816976803|gb|LAMW01000180.1|" 1055899, \
 "*gi|816977117|gb|LAMW01000075.1|" 1078413, \
 "gi|816977339|gb|LAMW01000001.1|" 1095351, \
 "*gi|816977294|gb|LAMW01000016.1|" 1111801, \
 "gi|816976815|gb|LAMW01000176.1|" 1122629, \
 "gi|816977024|gb|LAMW01000106.1|" 1149230, \
 "gi|816977009|gb|LAMW01000111.1|" 1162362, \
 "gi|816976964|gb|LAMW01000126.1|" 1173319, \
 "gi|816977318|gb|LAMW01000008.1|" 1186976, \
 "*gi|816976776|gb|LAMW01000189.1|" 1211658, \
 "gi|816977223|gb|LAMW01000040.1|" 1223079, \
 "gi|816977264|gb|LAMW01000026.1|" 1233234, \
 "gi|816976967|gb|LAMW01000125.1|" 1248013, \
 "gi|816976895|gb|LAMW01000149.1|" 1263903, \
 "gi|816976901|gb|LAMW01000147.1|" 1284860, \
 "*gi|816976818|gb|LAMW01000175.1|" 1311389, \
 "gi|816977247|gb|LAMW01000032.1|" 1326252, \
 "*gi|816976809|gb|LAMW01000178.1|" 1338941, \
 "*gi|816977027|gb|LAMW01000105.1|" 1352509, \
 "*gi|816977330|gb|LAMW01000004.1|" 1372162, \
 "*gi|816976851|gb|LAMW01000164.1|" 1383925, \
 "*gi|816977190|gb|LAMW01000051.1|" 1401666, \
 "gi|816977063|gb|LAMW01000093.1|" 1417143, \
 "*gi|816977090|gb|LAMW01000084.1|" 1427229, \
 "*gi|816977306|gb|LAMW01000012.1|" 1438644, \
 "*gi|816977018|gb|LAMW01000108.1|" 1449533, \
 "gi|816976931|gb|LAMW01000137.1|" 1468566, \
 "*gi|816977244|gb|LAMW01000033.1|" 1485890, \
 "*gi|816976922|gb|LAMW01000140.1|" 1499044, \
 "gi|816977051|gb|LAMW01000097.1|" 1512197, \
 "*gi|816976860|gb|LAMW01000161.1|" 1528421, \
 "*gi|816976797|gb|LAMW01000182.1|" 1542764, \
 "gi|816977146|gb|LAMW01000066.1|" 1564428, \
 "gi|816977006|gb|LAMW01000112.1|" 1577970, \
 "*gi|816977291|gb|LAMW01000017.1|" 1588768, \
 "gi|816977229|gb|LAMW01000038.1|" 1598956, \
 "*gi|816976883|gb|LAMW01000153.1|" 1611502, \
 "*gi|816977273|gb|LAMW01000023.1|" 1622684, \
 "gi|816976836|gb|LAMW01000169.1|" 1659409, \
 "gi|816976907|gb|LAMW01000145.1|" 1682987, \
 "*gi|816976782|gb|LAMW01000187.1|" 1696492, \
 "gi|816977282|gb|LAMW01000020.1|" 1708470, \
 "*gi|816976886|gb|LAMW01000152.1|" 1727978, \
 "*gi|816977066|gb|LAMW01000092.1|" 1743619, \
 "gi|816976952|gb|LAMW01000130.1|" 1762404, \
 "*gi|816976833|gb|LAMW01000170.1|" 1772727, \
 "gi|816977087|gb|LAMW01000085.1|" 1786313, \
 "*gi|816977124|gb|LAMW01000072.1|" 1797021, \
 "gi|816977134|gb|LAMW01000069.1|" 1815912, \
 "*gi|816977096|gb|LAMW01000082.1|" 1836520, \
 "gi|816976758|gb|LAMW01000195.1|" 1846955, \
 "gi|816977151|gb|LAMW01000064.1|" 1859404, \
 "*gi|816976892|gb|LAMW01000150.1|" 1882384, \
 "gi|816977170|gb|LAMW01000058.1|" 1893061, \
 "gi|816977315|gb|LAMW01000009.1|" 1911444, \
 "gi|816977021|gb|LAMW01000107.1|" 1924942, \
 "*gi|816976773|gb|LAMW01000190.1|" 1942162, \
 "*gi|816977211|gb|LAMW01000044.1|" 1958777, \
 "*gi|816976761|gb|LAMW01000194.1|" 1970072, \
 "gi|816976872|gb|LAMW01000157.1|" 1985187, \
 "*gi|816976994|gb|LAMW01000116.1|" 1996259, \
 "*gi|816977075|gb|LAMW01000089.1|" 2019697, \
 "gi|816977069|gb|LAMW01000091.1|" 2041007, \
 "gi|816977205|gb|LAMW01000046.1|" 2059469, \
 "gi|816977321|gb|LAMW01000007.1|" 2074755, \
 "gi|816977300|gb|LAMW01000014.1|" 2099263, \
 "*gi|816977048|gb|LAMW01000098.1|" 2117182, \
 "*gi|816976755|gb|LAMW01000196.1|" 2132498, \
 "*gi|816976800|gb|LAMW01000181.1|" 2149443, \
 "gi|816977297|gb|LAMW01000015.1|" 2165365, \
 "gi|816977093|gb|LAMW01000083.1|" 2181909, \
 "*gi|816977072|gb|LAMW01000090.1|" 2194463, \
 "gi|816977039|gb|LAMW01000101.1|" 2219480, \
 "*gi|816977193|gb|LAMW01000050.1|" 2229633, \
 "gi|816976940|gb|LAMW01000134.1|" 2252449, \
 "gi|816977060|gb|LAMW01000094.1|" 2262469, \
 "*gi|816976830|gb|LAMW01000171.1|" 2284186, \
 "*gi|816976946|gb|LAMW01000132.1|" 2298219, \
 "*gi|816977303|gb|LAMW01000013.1|" 2311815, \
 "*gi|816977214|gb|LAMW01000043.1|" 2332340, \
 "gi|816976740|gb|LAMW01000201.1|" 2342340, \
 "*gi|816977127|gb|LAMW01000071.1|" 2361663, \
 "*gi|816977181|gb|LAMW01000054.1|" 2374604, \
 "*gi|816977012|gb|LAMW01000110.1|" 2392756, \
 "gi|816977136|gb|LAMW01000068.1|" 2407032, \
 "*gi|816977184|gb|LAMW01000053.1|" 2435044, \
 "gi|816976925|gb|LAMW01000139.1|" 2450441, \
 "*gi|816976934|gb|LAMW01000136.1|" 2462629, \
 "*gi|816976937|gb|LAMW01000135.1|" 2476905, \
 "gi|816976949|gb|LAMW01000131.1|" 2490101, \
 "*gi|816976845|gb|LAMW01000166.1|" 2502277, \
 "*gi|816977324|gb|LAMW01000006.1|" 2517432, \
 "*gi|816977250|gb|LAMW01000031.1|" 2527455, \
 "gi|816976991|gb|LAMW01000117.1|" 2541055, \
 "*gi|816976943|gb|LAMW01000133.1|" 2569468, \
 "*gi|816977033|gb|LAMW01000103.1|" 2585637, \
 "gi|816977057|gb|LAMW01000095.1|" 2595691, \
 "gi|816977267|gb|LAMW01000025.1|" 2607490, \
 "*gi|816977036|gb|LAMW01000102.1|" 2617603, \
 "*gi|816976743|gb|LAMW01000200.1|" 2635579, \
 "*gi|816977154|gb|LAMW01000063.1|" 2645993, \
 "gi|816977327|gb|LAMW01000005.1|" 2668528, \
 "*gi|816976857|gb|LAMW01000162.1|" 2682528, \
 "*gi|816976839|gb|LAMW01000168.1|" 2694936, \
 "*gi|816977084|gb|LAMW01000086.1|" 2707410, \
 "gi|816977270|gb|LAMW01000024.1|" 2721693, \
 "gi|816976955|gb|LAMW01000129.1|" 2749964, \
 "*gi|816976749|gb|LAMW01000198.1|" 2760719, \
 "gi|816977309|gb|LAMW01000011.1|" 2773387, \
 "gi|816977187|gb|LAMW01000052.1|" 2786774, \
 "gi|816977081|gb|LAMW01000087.1|" 2815231, \
 "gi|816976889|gb|LAMW01000151.1|" 2839763, \
 "*gi|816977045|gb|LAMW01000099.1|" 2867558, \
 "gi|816976827|gb|LAMW01000172.1|" 2880476, \
 "gi|816976791|gb|LAMW01000184.1|" 2906126, \
 "gi|816977015|gb|LAMW01000109.1|" 2919027, \
 "*gi|816976854|gb|LAMW01000163.1|" 2931214, \
 "*gi|816976866|gb|LAMW01000159.1|" 2943634, \
 "gi|816977144|gb|LAMW01000067.1|" 2954564, \
 "gi|816976824|gb|LAMW01000173.1|" 2975655, \
 "gi|816977285|gb|LAMW01000019.1|" 2996265, \
 "*gi|816977102|gb|LAMW01000080.1|" 3012678, \
 "gi|816977078|gb|LAMW01000088.1|" 3024009, \
 "gi|816976764|gb|LAMW01000193.1|" 3039325, \
 "gi|816976982|gb|LAMW01000120.1|" 3054431, \
 "gi|816977000|gb|LAMW01000114.1|" 3068722, \
 "gi|816976973|gb|LAMW01000123.1|" 3082947, \
 "gi|816976770|gb|LAMW01000191.1|" 3096954, \
 "gi|816976985|gb|LAMW01000119.1|" 3109923, \
 "gi|816977114|gb|LAMW01000076.1|" 3122038, \
 "gi|816977262|gb|LAMW01000027.1|" 3133849, \
 "gi|816977199|gb|LAMW01000048.1|" 3145590, \
 "gi|816976970|gb|LAMW01000124.1|" 3156892, \
 "gi|816977108|gb|LAMW01000078.1|" 3168012, \
 "gi|816976806|gb|LAMW01000179.1|" 3179100, \
 "gi|816977226|gb|LAMW01000039.1|" 3189905, \
 "gi|816977279|gb|LAMW01000021.1|" 3200158, \
 "gi|816977259|gb|LAMW01000028.1|" 3213322, \
 "gi|816977161|gb|LAMW01000061.1|" 3230674, \
 "gi|816976919|gb|LAMW01000141.1|" 3246867, \
 "gi|816977336|gb|LAMW01000002.1|" 3259328, \
 "gi|816976794|gb|LAMW01000183.1|" 3274767, \
 "" 3286556 \
)
set ytics ( \
 "*JAAHHT010000106.1" 1, \
 "JAAHHT010000134.1" 30107, \
 "JAAHHT010000146.1" 55630, \
 "*JAAHHT010000309.1" 80272, \
 "JAAHHT010000160.1" 90074, \
 "JAAHHT010000075.1" 112438, \
 "*JAAHHT010000293.1" 150182, \
 "JAAHHT010000122.1" 160970, \
 "JAAHHT010000351.1" 187589, \
 "JAAHHT010000126.1" 194969, \
 "*JAAHHT010000117.1" 221195, \
 "JAAHHT010000138.1" 249052, \
 "*JAAHHT010000082.1" 274374, \
 "*JAAHHT010000034.1" 309537, \
 "*JAAHHT010000008.1" 373795, \
 "JAAHHT010000110.1" 480246, \
 "*JAAHHT010000027.1" 509173, \
 "JAAHHT010000400.1" 616109, \
 "*JAAHHT010000017.1" 621204, \
 "JAAHHT010000139.1" 821419, \
 "JAAHHT010000233.1" 846449, \
 "*JAAHHT010000028.1" 861834, \
 "*JAAHHT010000053.1" 963347, \
 "JAAHHT010000219.1" 1009095, \
 "*JAAHHT010000013.1" 1025419, \
 "JAAHHT010000254.1" 1111872, \
 "JAAHHT010000195.1" 1125085, \
 "*JAAHHT010000047.1" 1142934, \
 "JAAHHT010000010.1" 1190994, \
 "JAAHHT010000087.1" 1249223, \
 "JAAHHT010000172.1" 1283315, \
 "*JAAHHT010000016.1" 1304404, \
 "JAAHHT010000153.1" 1435999, \
 "JAAHHT010000057.1" 1458968, \
 "*JAAHHT010000045.1" 1503174, \
 "JAAHHT010000070.1" 1554538, \
 "JAAHHT010000068.1" 1594027, \
 "*JAAHHT010000249.1" 1633858, \
 "JAAHHT010000088.1" 1647443, \
 "*JAAHHT010000061.1" 1681168, \
 "*JAAHHT010000129.1" 1723562, \
 "JAAHHT010000097.1" 1749513, \
 "*JAAHHT010000189.1" 1781474, \
 "*JAAHHT010000125.1" 1800305, \
 "JAAHHT010000119.1" 1826704, \
 "*JAAHHT010000180.1" 1853866, \
 "*JAAHHT010000350.1" 1873932, \
 "JAAHHT010000022.1" 1881374, \
 "JAAHHT010000190.1" 1890782, \
 "*JAAHHT010000434.1" 1909495, \
 "*JAAHHT010000278.1" 1913197, \
 "*JAAHHT010000402.1" 1924664, \
 "JAAHHT010000186.1" 1929650, \
 "JAAHHT010000471.1" 1948894, \
 "JAAHHT010000015.1" 1951543, \
 "*JAAHHT010000214.1" 1987536, \
 "*JAAHHT010000223.1" 2004120, \
 "JAAHHT010000116.1" 2020202, \
 "JAAHHT010000200.1" 2048268, \
 "JAAHHT010000175.1" 2065714, \
 "*JAAHHT010000158.1" 2086351, \
 "JAAHHT010000339.1" 2108940, \
 "*JAAHHT010000379.1" 2116824, \
 "JAAHHT010000049.1" 2122807, \
 "*JAAHHT010000401.1" 2170527, \
 "JAAHHT010000221.1" 2175573, \
 "JAAHHT010000089.1" 2191770, \
 "JAAHHT010000271.1" 2225319, \
 "JAAHHT010000279.1" 2237300, \
 "*JAAHHT010000040.1" 2248759, \
 "JAAHHT010000210.1" 2308329, \
 "*JAAHHT010000035.1" 2325090, \
 "JAAHHT010000256.1" 2389188, \
 "JAAHHT010000043.1" 2402350, \
 "JAAHHT010000171.1" 2458067, \
 "*JAAHHT010000062.1" 2479243, \
 "*JAAHHT010000064.1" 2521094, \
 "JAAHHT010000182.1" 2562467, \
 "*JAAHHT010000019.1" 2581904, \
 "*JAAHHT010000124.1" 2710341, \
 "JAAHHT010000133.1" 2736800, \
 "JAAHHT010000202.1" 2762389, \
 "JAAHHT010000283.1" 2779687, \
 "*JAAHHT010000427.1" 2790886, \
 "JAAHHT010000224.1" 2794963, \
 "JAAHHT010000252.1" 2810726, \
 "JAAHHT010000154.1" 2824130, \
 "*JAAHHT010000318.1" 2847018, \
 "*JAAHHT010000162.1" 2856023, \
 "JAAHHT010000220.1" 2878104, \
 "JAAHHT010000397.1" 2894362, \
 "JAAHHT010000399.1" 2899657, \
 "JAAHHT010000054.1" 2904833, \
 "*JAAHHT010000216.1" 2950521, \
 "JAAHHT010000292.1" 2967088, \
 "JAAHHT010000163.1" 2977906, \
 "JAAHHT010000367.1" 2999978, \
 "JAAHHT010000020.1" 3006393, \
 "*JAAHHT010000038.1" 3043544, \
 "*JAAHHT010000031.1" 3105069, \
 "JAAHHT010000050.1" 3177224, \
 "JAAHHT010000007.1" 3224893, \
 "JAAHHT010000046.1" 3369850, \
 "*JAAHHT010000081.1" 3419328, \
 "JAAHHT010000002.1" 3454915, \
 "*JAAHHT010000181.1" 3523451, \
 "JAAHHT010000132.1" 3543017, \
 "JAAHHT010000033.1" 3568607, \
 "JAAHHT010000074.1" 3634843, \
 "*JAAHHT010000188.1" 3672771, \
 "JAAHHT010000229.1" 3691876, \
 "JAAHHT010000130.1" 3707450, \
 "*JAAHHT010000198.1" 3733379, \
 "JAAHHT010000204.1" 3750993, \
 "*JAAHHT010000251.1" 3768097, \
 "*JAAHHT010000001.1" 3781512, \
 "*JAAHHT010000018.1" 4016115, \
 "JAAHHT010000239.1" 4054194, \
 "JAAHHT010000288.1" 4068838, \
 "JAAHHT010000191.1" 4079718, \
 "JAAHHT010000264.1" 4098264, \
 "*JAAHHT010000073.1" 4110827, \
 "*JAAHHT010000205.1" 4149108, \
 "JAAHHT010000469.1" 4166072, \
 "JAAHHT010000011.1" 4168778, \
 "*JAAHHT010000072.1" 4240340, \
 "*JAAHHT010000004.1" 4279046, \
 "*JAAHHT010000036.1" 4359138, \
 "*JAAHHT010000083.1" 4421880, \
 "JAAHHT010000266.1" 4456834, \
 "JAAHHT010000437.1" 4469024, \
 "*JAAHHT010000128.1" 4472628, \
 "*JAAHHT010000490.1" 4498625, \
 "*JAAHHT010000131.1" 4500637, \
 "*JAAHHT010000085.1" 4526254, \
 "*JAAHHT010000482.1" 4560991, \
 "JAAHHT010000446.1" 4563176, \
 "*JAAHHT010000259.1" 4566422, \
 "*JAAHHT010000021.1" 4579416, \
 "*JAAHHT010000479.1" 4611567, \
 "*JAAHHT010000005.1" 4613816, \
 "JAAHHT010000300.1" 4699693, \
 "JAAHHT010000364.1" 4710026, \
 "JAAHHT010000137.1" 4716632, \
 "*JAAHHT010000317.1" 4741957, \
 "*JAAHHT010000065.1" 4751103, \
 "JAAHHT010000091.1" 4792441, \
 "JAAHHT010000164.1" 4825542, \
 "*JAAHHT010000023.1" 4847554, \
 "*JAAHHT010000103.1" 4876050, \
 "JAAHHT010000009.1" 4906673, \
 "JAAHHT010000140.1" 4985498, \
 "*JAAHHT010000238.1" 5010453, \
 "JAAHHT010000323.1" 5025169, \
 "*JAAHHT010000121.1" 5033891, \
 "*JAAHHT010000051.1" 5060710, \
 "*JAAHHT010000272.1" 5106942, \
 "JAAHHT010000060.1" 5118730, \
 "JAAHHT010000366.1" 5161869, \
 "*JAAHHT010000115.1" 5168314, \
 "*JAAHHT010000118.1" 5196625, \
 "JAAHHT010000340.1" 5224192, \
 "JAAHHT010000120.1" 5231991, \
 "JAAHHT010000167.1" 5258983, \
 "JAAHHT010000174.1" 5280520, \
 "JAAHHT010000148.1" 5301252, \
 "JAAHHT010000322.1" 5325356, \
 "JAAHHT010000307.1" 5334148, \
 "JAAHHT010000398.1" 5344239, \
 "JAAHHT010000384.1" 5349471, \
 "JAAHHT010000456.1" 5355290, \
 "*JAAHHT010000069.1" 5358374, \
 "*JAAHHT010000080.1" 5397900, \
 "JAAHHT010000326.1" 5433564, \
 "*JAAHHT010000184.1" 5442038, \
 "JAAHHT010000342.1" 5461363, \
 "JAAHHT010000369.1" 5469075, \
 "JAAHHT010000382.1" 5475422, \
 "JAAHHT010000165.1" 5481298, \
 "JAAHHT010000111.1" 5503293, \
 "*JAAHHT010000380.1" 5532118, \
 "*JAAHHT010000378.1" 5538078, \
 "JAAHHT010000459.1" 5544097, \
 "JAAHHT010000248.1" 5547064, \
 "*JAAHHT010000127.1" 5560770, \
 "*JAAHHT010000071.1" 5586943, \
 "JAAHHT010000149.1" 5625932, \
 "JAAHHT010000098.1" 5649309, \
 "JAAHHT010000063.1" 5680712, \
 "JAAHHT010000092.1" 5722411, \
 "JAAHHT010000056.1" 5755469, \
 "JAAHHT010000041.1" 5800336, \
 "JAAHHT010000094.1" 5859785, \
 "*JAAHHT010000394.1" 5892363, \
 "*JAAHHT010000463.1" 5897686, \
 "JAAHHT010000303.1" 5900584, \
 "JAAHHT010000241.1" 5910755, \
 "*JAAHHT010000301.1" 5924995, \
 "JAAHHT010000473.1" 5935250, \
 "JAAHHT010000168.1" 5937853, \
 "JAAHHT010000030.1" 5959111, \
 "JAAHHT010000178.1" 6031628, \
 "JAAHHT010000104.1" 6051780, \
 "*JAAHHT010000102.1" 6082167, \
 "*JAAHHT010000108.1" 6113018, \
 "*JAAHHT010000177.1" 6142761, \
 "JAAHHT010000012.1" 6163073, \
 "JAAHHT010000099.1" 6204192, \
 "*JAAHHT010000316.1" 6235473, \
 "*JAAHHT010000495.1" 6244809, \
 "JAAHHT010000329.1" 6246730, \
 "*JAAHHT010000282.1" 6255040, \
 "JAAHHT010000381.1" 6266324, \
 "*JAAHHT010000093.1" 6272265, \
 "*JAAHHT010000328.1" 6305061, \
 "*JAAHHT010000361.1" 6313437, \
 "JAAHHT010000055.1" 6320193, \
 "JAAHHT010000079.1" 6365780, \
 "JAAHHT010000144.1" 6401551, \
 "JAAHHT010000506.1" 6426336, \
 "*JAAHHT010000440.1" 6428084, \
 "JAAHHT010000147.1" 6431642, \
 "*JAAHHT010000026.1" 6455829, \
 "JAAHHT010000291.1" 6566825, \
 "*JAAHHT010000161.1" 6577676, \
 "*JAAHHT010000374.1" 6600029, \
 "*JAAHHT010000095.1" 6606218, \
 "*JAAHHT010000489.1" 6638723, \
 "*JAAHHT010000503.1" 6640743, \
 "*JAAHHT010000170.1" 6642535, \
 "*JAAHHT010000096.1" 6663720, \
 "*JAAHHT010000176.1" 6696068, \
 "JAAHHT010000370.1" 6716524, \
 "JAAHHT010000044.1" 6722867, \
 "*JAAHHT010000109.1" 6776452, \
 "JAAHHT010000359.1" 6805399, \
 "JAAHHT010000313.1" 6812241, \
 "*JAAHHT010000306.1" 6821713, \
 "JAAHHT010000435.1" 6831809, \
 "*JAAHHT010000208.1" 6835497, \
 "*JAAHHT010000123.1" 6852353, \
 "*JAAHHT010000406.1" 6878957, \
 "*JAAHHT010000487.1" 6883791, \
 "*JAAHHT010000531.1" 6885849, \
 "*JAAHHT010000247.1" 6887105, \
 "JAAHHT010000169.1" 6900814, \
 "*JAAHHT010000311.1" 6922034, \
 "*JAAHHT010000228.1" 6931731, \
 "JAAHHT010000245.1" 6947330, \
 "JAAHHT010000289.1" 6961195, \
 "*JAAHHT010000453.1" 6972057, \
 "*JAAHHT010000425.1" 6975208, \
 "*JAAHHT010000520.1" 6979347, \
 "*JAAHHT010000215.1" 6980818, \
 "*JAAHHT010000078.1" 6997388, \
 "JAAHHT010000166.1" 7033283, \
 "*JAAHHT010000207.1" 7055274, \
 "JAAHHT010000411.1" 7072200, \
 "JAAHHT010000039.1" 7076915, \
 "*JAAHHT010000357.1" 7137412, \
 "JAAHHT010000048.1" 7144347, \
 "JAAHHT010000273.1" 7192337, \
 "JAAHHT010000150.1" 7204083, \
 "*JAAHHT010000003.1" 7227403, \
 "*JAAHHT010000284.1" 7336449, \
 "JAAHHT010000356.1" 7347632, \
 "JAAHHT010000145.1" 7354636, \
 "JAAHHT010000305.1" 7379371, \
 "JAAHHT010000077.1" 7389488, \
 "JAAHHT010000525.1" 7426590, \
 "*JAAHHT010000114.1" 7427919, \
 "JAAHHT010000196.1" 7456293, \
 "*JAAHHT010000346.1" 7474035, \
 "JAAHHT010000025.1" 7481626, \
 "*JAAHHT010000212.1" 7607486, \
 "*JAAHHT010000173.1" 7624113, \
 "*JAAHHT010000006.1" 7644985, \
 "JAAHHT010000101.1" 7763677, \
 "*JAAHHT010000143.1" 7794635, \
 "JAAHHT010000076.1" 7819493, \
 "JAAHHT010000067.1" 7857167, \
 "*JAAHHT010000476.1" 7897500, \
 "JAAHHT010000410.1" 7900060, \
 "JAAHHT010000227.1" 7904786, \
 "JAAHHT010000058.1" 7920418, \
 "*JAAHHT010000090.1" 7964604, \
 "JAAHHT010000206.1" 7997708, \
 "JAAHHT010000260.1" 8014637, \
 "JAAHHT010000348.1" 8027552, \
 "JAAHHT010000516.1" 8035100, \
 "*JAAHHT010000280.1" 8036614, \
 "*JAAHHT010000324.1" 8047991, \
 "JAAHHT010000084.1" 8056595, \
 "*JAAHHT010000042.1" 8091501, \
 "*JAAHHT010000014.1" 8149002, \
 "JAAHHT010000203.1" 8207283, \
 "JAAHHT010000347.1" 8224459, \
 "*JAAHHT010000157.1" 8232012, \
 "JAAHHT010000561.1" 8254818, \
 "JAAHHT010000500.1" 8255821, \
 "JAAHHT010000553.1" 8257672, \
 "JAAHHT010000545.1" 8258709, \
 "JAAHHT010000066.1" 8259828, \
 "JAAHHT010000032.1" 8300906, \
 "JAAHHT010000353.1" 8372603, \
 "JAAHHT010000257.1" 8379821, \
 "JAAHHT010000389.1" 8392977, \
 "JAAHHT010000352.1" 8398522, \
 "JAAHHT010000308.1" 8405878, \
 "JAAHHT010000396.1" 8415818, \
 "JAAHHT010000496.1" 8421141, \
 "JAAHHT010000385.1" 8423034, \
 "JAAHHT010000310.1" 8428800, \
 "JAAHHT010000466.1" 8438557, \
 "JAAHHT010000275.1" 8441405, \
 "JAAHHT010000423.1" 8453024, \
 "JAAHHT010000187.1" 8457278, \
 "JAAHHT010000334.1" 8476421, \
 "JAAHHT010000037.1" 8484514, \
 "JAAHHT010000304.1" 8546925, \
 "JAAHHT010000421.1" 8557071, \
 "JAAHHT010000246.1" 8561359, \
 "JAAHHT010000193.1" 8575224, \
 "JAAHHT010000559.1" 8593507, \
 "JAAHHT010000515.1" 8594515, \
 "JAAHHT010000029.1" 8596068, \
 "JAAHHT010000414.1" 8682783, \
 "JAAHHT010000365.1" 8687439, \
 "JAAHHT010000136.1" 8693954, \
 "JAAHHT010000331.1" 8719379, \
 "JAAHHT010000226.1" 8727528, \
 "JAAHHT010000321.1" 8743199, \
 "JAAHHT010000253.1" 8752066, \
 "JAAHHT010000314.1" 8765304, \
 "JAAHHT010000554.1" 8774747, \
 "JAAHHT010000338.1" 8775783, \
 "JAAHHT010000112.1" 8783689, \
 "JAAHHT010000526.1" 8812406, \
 "JAAHHT010000231.1" 8813727, \
 "JAAHHT010000192.1" 8829265, \
 "JAAHHT010000507.1" 8847714, \
 "JAAHHT010000416.1" 8849450, \
 "JAAHHT010000312.1" 8853982, \
 "JAAHHT010000325.1" 8863496, \
 "JAAHHT010000483.1" 8872042, \
 "JAAHHT010000556.1" 8874194, \
 "JAAHHT010000105.1" 8875215, \
 "JAAHHT010000426.1" 8905503, \
 "JAAHHT010000420.1" 8909616, \
 "JAAHHT010000290.1" 8913953, \
 "JAAHHT010000151.1" 8924812, \
 "JAAHHT010000261.1" 8948106, \
 "JAAHHT010000360.1" 8960899, \
 "JAAHHT010000470.1" 8967656, \
 "JAAHHT010000240.1" 8970307, \
 "JAAHHT010000344.1" 8984784, \
 "JAAHHT010000409.1" 8992407, \
 "JAAHHT010000462.1" 8997180, \
 "JAAHHT010000413.1" 9000094, \
 "JAAHHT010000530.1" 9004787, \
 "JAAHHT010000281.1" 9006060, \
 "JAAHHT010000521.1" 9017376, \
 "JAAHHT010000250.1" 9018797, \
 "JAAHHT010000475.1" 9032260, \
 "JAAHHT010000390.1" 9034835, \
 "JAAHHT010000464.1" 9040363, \
 "JAAHHT010000529.1" 9043233, \
 "JAAHHT010000493.1" 9044524, \
 "JAAHHT010000412.1" 9046475, \
 "JAAHHT010000422.1" 9051173, \
 "JAAHHT010000465.1" 9055429, \
 "JAAHHT010000468.1" 9058284, \
 "JAAHHT010000262.1" 9061030, \
 "JAAHHT010000267.1" 9073799, \
 "JAAHHT010000299.1" 9085994, \
 "JAAHHT010000546.1" 9096351, \
 "JAAHHT010000373.1" 9097470, \
 "JAAHHT010000296.1" 9103705, \
 "JAAHHT010000540.1" 9114249, \
 "JAAHHT010000194.1" 9115424, \
 "JAAHHT010000059.1" 9133643, \
 "JAAHHT010000543.1" 9177705, \
 "JAAHHT010000444.1" 9178855, \
 "JAAHHT010000519.1" 9182284, \
 "JAAHHT010000159.1" 9183758, \
 "JAAHHT010000455.1" 9206275, \
 "JAAHHT010000244.1" 9209385, \
 "JAAHHT010000086.1" 9223404, \
 "JAAHHT010000498.1" 9258079, \
 "JAAHHT010000535.1" 9259962, \
 "JAAHHT010000343.1" 9261208, \
 "JAAHHT010000499.1" 9268878, \
 "JAAHHT010000209.1" 9270731, \
 "JAAHHT010000319.1" 9287571, \
 "JAAHHT010000211.1" 9296554, \
 "JAAHHT010000432.1" 9313264, \
 "JAAHHT010000395.1" 9317049, \
 "JAAHHT010000295.1" 9322372, \
 "JAAHHT010000375.1" 9333037, \
 "JAAHHT010000445.1" 9339165, \
 "JAAHHT010000024.1" 9342551, \
 "JAAHHT010000358.1" 9346739, \
 "JAAHHT010000480.1" 9353603, \
 "JAAHHT010000234.1" 9355824, \
 "JAAHHT010000532.1" 9370903, \
 "JAAHHT010000222.1" 9372158, \
 "JAAHHT010000428.1" 9388274, \
 "JAAHHT010000388.1" 9392337, \
 "JAAHHT010000492.1" 9397896, \
 "JAAHHT010000512.1" 9399855, \
 "JAAHHT010000407.1" 9401444, \
 "JAAHHT010000477.1" 9406213, \
 "JAAHHT010000458.1" 9408750, \
 "JAAHHT010000107.1" 9411717, \
 "JAAHHT010000518.1" 9441652, \
 "JAAHHT010000335.1" 9443133, \
 "JAAHHT010000502.1" 9451155, \
 "JAAHHT010000418.1" 9452950, \
 "JAAHHT010000135.1" 9457374, \
 "JAAHHT010000441.1" 9482800, \
 "JAAHHT010000302.1" 9486345, \
 "JAAHHT010000243.1" 9496555, \
 "JAAHHT010000217.1" 9510688, \
 "JAAHHT010000501.1" 9527062, \
 "JAAHHT010000443.1" 9528904, \
 "JAAHHT010000294.1" 9532365, \
 "JAAHHT010000392.1" 9543142, \
 "JAAHHT010000442.1" 9548633, \
 "JAAHHT010000509.1" 9552116, \
 "JAAHHT010000333.1" 9553786, \
 "JAAHHT010000393.1" 9561888, \
 "JAAHHT010000486.1" 9567259, \
 "JAAHHT010000514.1" 9569335, \
 "JAAHHT010000457.1" 9570907, \
 "JAAHHT010000542.1" 9573889, \
 "JAAHHT010000258.1" 9575040, \
 "JAAHHT010000419.1" 9588142, \
 "JAAHHT010000404.1" 9592482, \
 "JAAHHT010000332.1" 9597374, \
 "JAAHHT010000315.1" 9605483, \
 "JAAHHT010000539.1" 9614851, \
 "JAAHHT010000179.1" 9616038, \
 "JAAHHT010000417.1" 9636118, \
 "JAAHHT010000485.1" 9640629, \
 "JAAHHT010000341.1" 9642756, \
 "JAAHHT010000230.1" 9650492, \
 "JAAHHT010000562.1" 9666043, \
 "JAAHHT010000274.1" 9667045, \
 "JAAHHT010000550.1" 9678732, \
 "JAAHHT010000268.1" 9679824, \
 "JAAHHT010000474.1" 9691960, \
 "JAAHHT010000362.1" 9694562, \
 "JAAHHT010000439.1" 9701263, \
 "JAAHHT010000213.1" 9704833, \
 "JAAHHT010000287.1" 9721432, \
 "JAAHHT010000488.1" 9732325, \
 "JAAHHT010000549.1" 9734364, \
 "JAAHHT010000270.1" 9735457, \
 "JAAHHT010000505.1" 9747477, \
 "JAAHHT010000327.1" 9749243, \
 "JAAHHT010000430.1" 9757626, \
 "JAAHHT010000298.1" 9761506, \
 "JAAHHT010000336.1" 9771902, \
 "JAAHHT010000448.1" 9779900, \
 "JAAHHT010000447.1" 9783100, \
 "JAAHHT010000541.1" 9786302, \
 "JAAHHT010000491.1" 9787458, \
 "JAAHHT010000478.1" 9789462, \
 "JAAHHT010000255.1" 9791860, \
 "JAAHHT010000377.1" 9805056, \
 "JAAHHT010000544.1" 9811080, \
 "JAAHHT010000405.1" 9812205, \
 "JAAHHT010000538.1" 9817039, \
 "JAAHHT010000523.1" 9818231, \
 "JAAHHT010000534.1" 9819612, \
 "JAAHHT010000533.1" 9820862, \
 "JAAHHT010000286.1" 9822114, \
 "JAAHHT010000235.1" 9833067, \
 "JAAHHT010000386.1" 9848145, \
 "JAAHHT010000265.1" 9853735, \
 "JAAHHT010000355.1" 9866240, \
 "JAAHHT010000547.1" 9873294, \
 "JAAHHT010000201.1" 9874412, \
 "JAAHHT010000403.1" 9891840, \
 "JAAHHT010000276.1" 9896783, \
 "JAAHHT010000263.1" 9908342, \
 "JAAHHT010000433.1" 9921094, \
 "JAAHHT010000142.1" 9924816, \
 "JAAHHT010000383.1" 9949749, \
 "JAAHHT010000504.1" 9955618, \
 "JAAHHT010000337.1" 9957401, \
 "JAAHHT010000242.1" 9965352, \
 "JAAHHT010000155.1" 9979523, \
 "JAAHHT010000052.1" 10002373, \
 "JAAHHT010000451.1" 10048489, \
 "JAAHHT010000552.1" 10051659, \
 "JAAHHT010000548.1" 10052730, \
 "JAAHHT010000354.1" 10053824, \
 "JAAHHT010000454.1" 10060984, \
 "JAAHHT010000376.1" 10064115, \
 "JAAHHT010000481.1" 10070211, \
 "JAAHHT010000199.1" 10072429, \
 "JAAHHT010000141.1" 10089914, \
 "JAAHHT010000436.1" 10114853, \
 "JAAHHT010000494.1" 10118490, \
 "JAAHHT010000555.1" 10120439, \
 "JAAHHT010000497.1" 10121468, \
 "JAAHHT010000372.1" 10123358, \
 "JAAHHT010000387.1" 10129611, \
 "JAAHHT010000510.1" 10135194, \
 "JAAHHT010000345.1" 10136806, \
 "JAAHHT010000156.1" 10144402, \
 "JAAHHT010000438.1" 10167237, \
 "JAAHHT010000527.1" 10170819, \
 "JAAHHT010000517.1" 10172134, \
 "JAAHHT010000368.1" 10173644, \
 "JAAHHT010000460.1" 10180014, \
 "JAAHHT010000537.1" 10182965, \
 "JAAHHT010000524.1" 10184174, \
 "JAAHHT010000472.1" 10185548, \
 "JAAHHT010000452.1" 10188164, \
 "JAAHHT010000237.1" 10191326, \
 "JAAHHT010000285.1" 10206060, \
 "JAAHHT010000152.1" 10217233, \
 "JAAHHT010000536.1" 10240392, \
 "JAAHHT010000320.1" 10241621, \
 "JAAHHT010000297.1" 10250519, \
 "JAAHHT010000218.1" 10260994, \
 "JAAHHT010000513.1" 10277367, \
 "JAAHHT010000236.1" 10278944, \
 "JAAHHT010000232.1" 10293884, \
 "JAAHHT010000431.1" 10309357, \
 "JAAHHT010000185.1" 10313152, \
 "JAAHHT010000467.1" 10332431, \
 "JAAHHT010000508.1" 10335226, \
 "JAAHHT010000511.1" 10336898, \
 "JAAHHT010000371.1" 10338506, \
 "JAAHHT010000560.1" 10344821, \
 "JAAHHT010000450.1" 10345825, \
 "JAAHHT010000113.1" 10349010, \
 "JAAHHT010000484.1" 10377403, \
 "JAAHHT010000424.1" 10379531, \
 "JAAHHT010000558.1" 10383784, \
 "JAAHHT010000449.1" 10384796, \
 "JAAHHT010000363.1" 10387993, \
 "JAAHHT010000461.1" 10394668, \
 "JAAHHT010000408.1" 10397597, \
 "JAAHHT010000415.1" 10402396, \
 "JAAHHT010000183.1" 10407006, \
 "JAAHHT010000277.1" 10426362, \
 "JAAHHT010000330.1" 10437891, \
 "JAAHHT010000269.1" 10446123, \
 "JAAHHT010000197.1" 10458255, \
 "JAAHHT010000225.1" 10475977, \
 "JAAHHT010000522.1" 10491697, \
 "JAAHHT010000349.1" 10493100, \
 "JAAHHT010000557.1" 10500592, \
 "JAAHHT010000100.1" 10501609, \
 "JAAHHT010000528.1" 10532773, \
 "JAAHHT010000551.1" 10534080, \
 "JAAHHT010000429.1" 10535170, \
 "JAAHHT010000391.1" 10539225, \
 "" 10545295 \
)
set size 1,1
set grid
unset key
set border 0
set tics scale 0
set xlabel "REF"
set ylabel "QRY"
set format "%.0f"
set mouse format "%.0f"
set mouse mouseformat "[%.0f, %.0f]"
if(GPVAL_VERSION < 5) set mouse clipboardformat "[%.0f, %.0f]"
set xrange [1:3286556]
set yrange [1:10545295]
set style line 1  lt 1 lw 2 pt 6 ps 0.5
set style line 2  lt 3 lw 2 pt 6 ps 0.5
set style line 3  lt 2 lw 2 pt 6 ps 0.5
plot \
 "LAMW01_vs_JAAHHT01.fplot" title "FWD" w lp ls 1, \
 "LAMW01_vs_JAAHHT01.rplot" title "REV" w lp ls 2
