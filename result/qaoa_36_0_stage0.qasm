OPENQASM 2.0;
include "qelib1.inc";
qreg q[127];
creg c[127];
rzz(pi/4) q[106], q[107];
rzz(pi/4) q[82], q[83];
rzz(pi/4) q[121], q[120];
rzz(pi/4) q[126], q[125];
rzz(pi/4) q[122], q[123];
rzz(pi/4) q[101], q[102];
swap q[92], q[102];
swap q[100], q[101];
swap q[81], q[82];
swap q[107], q[108];
swap q[93], q[106];
swap q[112], q[126];
swap q[124], q[125];
swap q[121], q[122];
rzz(pi/4) q[81], q[80];
swap q[83], q[92];
swap q[100], q[110];
swap q[87], q[93];
swap q[108], q[112];
swap q[123], q[124];
rzz(pi/4) q[87], q[88];
rzz(pi/4) q[100], q[101];
rzz(pi/4) q[107], q[108];
swap q[110], q[118];
swap q[92], q[102];
swap q[112], q[126];
swap q[79], q[80];
swap q[102], q[103];
swap q[118], q[119];
swap q[86], q[87];
swap q[106], q[107];
swap q[125], q[126];
swap q[100], q[110];
swap q[79], q[91];
rzz(pi/4) q[100], q[99];
swap q[103], q[104];
swap q[85], q[86];
swap q[93], q[106];
swap q[110], q[118];
swap q[101], q[102];
swap q[91], q[98];
swap q[119], q[120];
rzz(pi/4) q[79], q[91];
swap q[104], q[111];
swap q[84], q[85];
swap q[87], q[93];
swap q[118], q[119];
swap q[100], q[110];
swap q[120], q[121];
rzz(pi/4) q[112], q[111];
rzz(pi/4) q[93], q[106];
rzz(pi/4) q[104], q[103];
rzz(pi/4) q[86], q[85];
rzz(pi/4) q[120], q[121];
rzz(pi/4) q[79], q[80];
rzz(pi/4) q[98], q[91];
swap q[83], q[84];
swap q[110], q[118];
rzz(pi/4) q[98], q[99];
rzz(pi/4) q[120], q[119];
rzz(pi/4) q[122], q[121];
swap q[104], q[105];
swap q[102], q[103];
swap q[86], q[87];
swap q[112], q[126];
swap q[78], q[79];
rzz(pi/4) q[122], q[123];
rzz(pi/4) q[105], q[106];
rzz(pi/4) q[102], q[101];
rzz(pi/4) q[126], q[112];
rzz(pi/4) q[78], q[79];
swap q[87], q[93];
rzz(pi/4) q[126], q[125];
rzz(pi/4) q[124], q[123];
rzz(pi/4) q[86], q[87];
rzz(pi/4) q[91], q[79];
swap q[100], q[101];
swap q[92], q[102];
swap q[106], q[107];
rzz(pi/4) q[125], q[124];
rzz(pi/4) q[101], q[102];
rzz(pi/4) q[87], q[88];
rzz(pi/4) q[93], q[106];
rzz(pi/4) q[86], q[85];
rzz(pi/4) q[80], q[79];
swap q[100], q[110];
rzz(pi/4) q[81], q[80];
rzz(pi/4) q[85], q[84];
rzz(pi/4) q[110], q[118];
swap q[100], q[101];
swap q[87], q[93];
rzz(pi/4) q[83], q[84];
rzz(pi/4) q[100], q[110];
rzz(pi/4) q[101], q[102];
rzz(pi/4) q[119], q[118];
rzz(pi/4) q[87], q[88];
rzz(pi/4) q[100], q[99];
swap q[102], q[103];
rzz(pi/4) q[92], q[102];
rzz(pi/4) q[104], q[103];
rzz(pi/4) q[102], q[101];
swap q[104], q[105];
rzz(pi/4) q[106], q[105];
swap q[104], q[111];
rzz(pi/4) q[104], q[105];
rzz(pi/4) q[111], q[112];
rzz(pi/4) q[106], q[107];
rzz(pi/4) q[108], q[112];

// measurement
measure q[120]->c[0];
measure q[86]->c[1];
measure q[87]->c[2];
measure q[93]->c[3];
measure q[119]->c[4];
measure q[122]->c[5];
measure q[81]->c[6];
measure q[108]->c[7];
measure q[111]->c[8];
measure q[106]->c[9];
measure q[100]->c[10];
measure q[85]->c[11];
measure q[92]->c[12];
measure q[83]->c[13];
measure q[78]->c[14];
measure q[110]->c[15];
measure q[126]->c[16];
measure q[121]->c[17];
measure q[107]->c[18];
measure q[98]->c[19];
measure q[80]->c[20];
measure q[102]->c[21];
measure q[101]->c[22];
measure q[91]->c[23];
measure q[104]->c[24];
measure q[88]->c[25];
measure q[118]->c[26];
measure q[125]->c[27];
measure q[124]->c[28];
measure q[105]->c[29];
measure q[79]->c[30];
measure q[112]->c[31];
measure q[123]->c[32];
measure q[103]->c[33];
measure q[99]->c[34];
measure q[84]->c[35];
