OPENQASM 2.0;
include "qelib1.inc";
qreg q[144];
creg c[144];
rz(pi/2) q[12];
sx q[12];
rz(pi/2) q[12];
cx q[12], q[24];
cx q[24], q[36];
cx q[36], q[48];
cx q[48], q[60];
cx q[60], q[72];
cx q[72], q[84];
cx q[84], q[85];
cx q[85], q[86];
cx q[86], q[98];
cx q[98], q[97];
cx q[97], q[96];
cx q[96], q[108];
cx q[108], q[120];
cx q[120], q[121];
cx q[121], q[109];
cx q[109], q[110];
cx q[110], q[111];
cx q[111], q[112];
cx q[112], q[124];
cx q[124], q[123];
cx q[123], q[122];
cx q[122], q[134];
cx q[134], q[135];
cx q[135], q[136];
cx q[136], q[137];
cx q[137], q[138];
cx q[138], q[139];
cx q[139], q[140];
cx q[140], q[141];
cx q[141], q[142];
cx q[142], q[130];
cx q[130], q[118];
cx q[118], q[106];
cx q[106], q[105];
cx q[105], q[93];
cx q[93], q[94];
cx q[94], q[95];
cx q[95], q[83];
cx q[83], q[71];
cx q[71], q[59];
cx q[59], q[58];
cx q[58], q[57];
cx q[57], q[45];
cx q[45], q[46];
cx q[46], q[34];
cx q[34], q[33];
cx q[33], q[21];
cx q[21], q[22];
cx q[22], q[10];
cx q[10], q[9];
cx q[9], q[8];
cx q[8], q[7];
cx q[7], q[6];
cx q[6], q[5];
cx q[5], q[4];
cx q[4], q[3];
cx q[3], q[2];
cx q[2], q[1];
cx q[1], q[13];
cx q[13], q[25];
cx q[25], q[37];
cx q[37], q[38];
cx q[38], q[39];
cx q[39], q[51];
cx q[51], q[50];
cx q[50], q[49];
cx q[49], q[61];
cx q[61], q[73];
cx q[73], q[74];
cx q[74], q[62];
cx q[62], q[63];
cx q[63], q[75];
cx q[75], q[87];
cx q[87], q[99];
cx q[99], q[100];
cx q[100], q[88];
cx q[88], q[76];
cx q[76], q[64];
cx q[64], q[52];
cx q[52], q[40];
cx q[40], q[41];
cx q[41], q[42];
cx q[42], q[54];
cx q[54], q[53];
cx q[53], q[65];
cx q[65], q[66];
cx q[66], q[78];
cx q[78], q[77];
cx q[77], q[89];
cx q[89], q[90];
cx q[90], q[102];
cx q[102], q[114];
cx q[114], q[113];
cx q[113], q[125];
cx q[125], q[126];
cx q[126], q[127];
cx q[127], q[128];
cx q[128], q[116];
cx q[116], q[115];
cx q[115], q[103];
cx q[103], q[104];
cx q[104], q[92];
cx q[92], q[80];
cx q[80], q[81];
cx q[81], q[82];
cx q[82], q[70];
cx q[70], q[69];
cx q[69], q[68];
cx q[68], q[56];
cx q[56], q[44];
cx q[44], q[32];
cx q[32], q[20];
cx q[20], q[19];
cx q[19], q[18];
cx q[18], q[17];
cx q[17], q[16];
cx q[16], q[15];
cx q[15], q[14];
cx q[14], q[26];
cx q[26], q[27];
cx q[27], q[28];
cx q[28], q[29];
cx q[29], q[30];
cx q[30], q[31];
cx q[31], q[43];
cx q[43], q[55];
cx q[55], q[67];
cx q[67], q[79];
cx q[79], q[91];

// measurement
measure q[12]->c[0];
measure q[24]->c[1];
measure q[36]->c[2];
measure q[48]->c[3];
measure q[60]->c[4];
measure q[72]->c[5];
measure q[84]->c[6];
measure q[85]->c[7];
measure q[86]->c[8];
measure q[98]->c[9];
measure q[97]->c[10];
measure q[96]->c[11];
measure q[108]->c[12];
measure q[120]->c[13];
measure q[121]->c[14];
measure q[109]->c[15];
measure q[110]->c[16];
measure q[111]->c[17];
measure q[112]->c[18];
measure q[124]->c[19];
measure q[123]->c[20];
measure q[122]->c[21];
measure q[134]->c[22];
measure q[135]->c[23];
measure q[136]->c[24];
measure q[137]->c[25];
measure q[138]->c[26];
measure q[139]->c[27];
measure q[140]->c[28];
measure q[141]->c[29];
measure q[142]->c[30];
measure q[130]->c[31];
measure q[118]->c[32];
measure q[106]->c[33];
measure q[105]->c[34];
measure q[93]->c[35];
measure q[94]->c[36];
measure q[95]->c[37];
measure q[83]->c[38];
measure q[71]->c[39];
measure q[59]->c[40];
measure q[58]->c[41];
measure q[57]->c[42];
measure q[45]->c[43];
measure q[46]->c[44];
measure q[34]->c[45];
measure q[33]->c[46];
measure q[21]->c[47];
measure q[22]->c[48];
measure q[10]->c[49];
measure q[9]->c[50];
measure q[8]->c[51];
measure q[7]->c[52];
measure q[6]->c[53];
measure q[5]->c[54];
measure q[4]->c[55];
measure q[3]->c[56];
measure q[2]->c[57];
measure q[1]->c[58];
measure q[13]->c[59];
measure q[25]->c[60];
measure q[37]->c[61];
measure q[38]->c[62];
measure q[39]->c[63];
measure q[51]->c[64];
measure q[50]->c[65];
measure q[49]->c[66];
measure q[61]->c[67];
measure q[73]->c[68];
measure q[74]->c[69];
measure q[62]->c[70];
measure q[63]->c[71];
measure q[75]->c[72];
measure q[87]->c[73];
measure q[99]->c[74];
measure q[100]->c[75];
measure q[88]->c[76];
measure q[76]->c[77];
measure q[64]->c[78];
measure q[52]->c[79];
measure q[40]->c[80];
measure q[41]->c[81];
measure q[42]->c[82];
measure q[54]->c[83];
measure q[53]->c[84];
measure q[65]->c[85];
measure q[66]->c[86];
measure q[78]->c[87];
measure q[77]->c[88];
measure q[89]->c[89];
measure q[90]->c[90];
measure q[102]->c[91];
measure q[114]->c[92];
measure q[113]->c[93];
measure q[125]->c[94];
measure q[126]->c[95];
measure q[127]->c[96];
measure q[128]->c[97];
measure q[116]->c[98];
measure q[115]->c[99];
measure q[103]->c[100];
measure q[104]->c[101];
measure q[92]->c[102];
measure q[80]->c[103];
measure q[81]->c[104];
measure q[82]->c[105];
measure q[70]->c[106];
measure q[69]->c[107];
measure q[68]->c[108];
measure q[56]->c[109];
measure q[44]->c[110];
measure q[32]->c[111];
measure q[20]->c[112];
measure q[19]->c[113];
measure q[18]->c[114];
measure q[17]->c[115];
measure q[16]->c[116];
measure q[15]->c[117];
measure q[14]->c[118];
measure q[26]->c[119];
measure q[27]->c[120];
measure q[28]->c[121];
measure q[29]->c[122];
measure q[30]->c[123];
measure q[31]->c[124];
measure q[43]->c[125];
measure q[55]->c[126];
measure q[67]->c[127];
measure q[79]->c[128];
measure q[91]->c[129];