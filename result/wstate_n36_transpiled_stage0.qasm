OPENQASM 2.0;
include "qelib1.inc";
qreg q[127];
creg c[127];
sx q[70];
sx q[74];
sx q[89];
sx q[88];
sx q[87];
sx q[93];
sx q[106];
sx q[107];
sx q[108];
sx q[112];
sx q[126];
sx q[125];
sx q[124];
sx q[123];
sx q[122];
sx q[121];
sx q[120];
sx q[119];
sx q[118];
sx q[110];
sx q[100];
sx q[101];
sx q[102];
sx q[92];
sx q[83];
sx q[84];
sx q[85];
sx q[73];
sx q[66];
sx q[65];
sx q[64];
sx q[63];
sx q[62];
sx q[61];
sx q[60];
x q[53];
rz(pi/4) q[70];
rz(0.61547971) q[74];
rz(pi/6) q[89];
rz(0.46364763) q[88];
rz(0.42053433) q[87];
rz(0.38759673) q[93];
rz(0.36136713) q[106];
rz(0.33983693) q[107];
rz(0.32175053) q[108];
rz(0.30627733) q[112];
rz(0.29284273) q[126];
rz(0.28103493) q[125];
rz(0.27054973) q[124];
rz(0.26115743) q[123];
rz(0.25268023) q[122];
rz(0.24497863) q[121];
rz(0.23794113) q[120];
rz(0.23147733) q[119];
rz(0.22551343) q[118];
rz(0.21998803) q[110];
rz(0.21484983) q[100];
rz(0.21005573) q[101];
rz(0.20556893) q[102];
rz(0.20135793) q[92];
rz(0.19739553) q[83];
rz(0.19365833) q[84];
rz(0.19012563) q[85];
rz(0.18677943) q[73];
rz(0.18360403) q[66];
rz(0.18058523) q[65];
rz(0.17771063) q[64];
rz(0.17496903) q[63];
rz(0.17235063) q[62];
rz(0.16984633) q[61];
rz(0.16744813) q[60];
sx q[70];
sx q[74];
sx q[89];
sx q[88];
sx q[87];
sx q[93];
sx q[106];
sx q[107];
sx q[108];
sx q[112];
sx q[126];
sx q[125];
sx q[124];
sx q[123];
sx q[122];
sx q[121];
sx q[120];
sx q[119];
sx q[118];
sx q[110];
sx q[100];
sx q[101];
sx q[102];
sx q[92];
sx q[83];
sx q[84];
sx q[85];
sx q[73];
sx q[66];
sx q[65];
sx q[64];
sx q[63];
sx q[62];
sx q[61];
sx q[60];
cx q[60], q[53];
sx q[60];
rz(0.16744813) q[60];
sx q[60];
cx q[61], q[60];
sx q[61];
cx q[60], q[53];
rz(0.16984633) q[61];
sx q[61];
cx q[62], q[61];
sx q[62];
cx q[61], q[60];
rz(0.17235063) q[62];
sx q[62];
cx q[63], q[62];
sx q[63];
cx q[62], q[61];
rz(0.17496903) q[63];
sx q[63];
cx q[64], q[63];
sx q[64];
cx q[63], q[62];
rz(0.17771063) q[64];
sx q[64];
cx q[65], q[64];
sx q[65];
cx q[64], q[63];
rz(0.18058523) q[65];
sx q[65];
cx q[66], q[65];
sx q[66];
cx q[65], q[64];
rz(0.18360403) q[66];
sx q[66];
cx q[73], q[66];
sx q[73];
cx q[66], q[65];
rz(0.18677943) q[73];
sx q[73];
cx q[85], q[73];
sx q[85];
cx q[73], q[66];
rz(0.19012563) q[85];
sx q[85];
cx q[84], q[85];
sx q[84];
cx q[85], q[73];
rz(0.19365833) q[84];
sx q[84];
cx q[83], q[84];
sx q[83];
cx q[84], q[85];
rz(0.19739553) q[83];
sx q[83];
cx q[92], q[83];
sx q[92];
cx q[83], q[84];
rz(0.20135793) q[92];
sx q[92];
cx q[102], q[92];
sx q[102];
cx q[92], q[83];
rz(0.20556893) q[102];
sx q[102];
cx q[101], q[102];
sx q[101];
cx q[102], q[92];
rz(0.21005573) q[101];
sx q[101];
cx q[100], q[101];
sx q[100];
cx q[101], q[102];
rz(0.21484983) q[100];
sx q[100];
cx q[110], q[100];
sx q[110];
cx q[100], q[101];
rz(0.21998803) q[110];
sx q[110];
cx q[118], q[110];
sx q[118];
cx q[110], q[100];
rz(0.22551343) q[118];
sx q[118];
cx q[119], q[118];
sx q[119];
cx q[118], q[110];
rz(0.23147733) q[119];
sx q[119];
cx q[120], q[119];
sx q[120];
cx q[119], q[118];
rz(0.23794113) q[120];
sx q[120];
cx q[121], q[120];
sx q[121];
cx q[120], q[119];
rz(0.24497863) q[121];
sx q[121];
cx q[122], q[121];
sx q[122];
cx q[121], q[120];
rz(0.25268023) q[122];
sx q[122];
cx q[123], q[122];
sx q[123];
cx q[122], q[121];
rz(0.26115743) q[123];
sx q[123];
cx q[124], q[123];
sx q[124];
cx q[123], q[122];
rz(0.27054973) q[124];
sx q[124];
cx q[125], q[124];
sx q[125];
cx q[124], q[123];
rz(0.28103493) q[125];
sx q[125];
cx q[126], q[125];
sx q[126];
cx q[125], q[124];
rz(0.29284273) q[126];
sx q[126];
cx q[112], q[126];
cx q[126], q[125];
sx q[112];
rz(0.30627733) q[112];
sx q[112];
cx q[108], q[112];
sx q[108];
cx q[112], q[126];
rz(0.32175053) q[108];
sx q[108];
cx q[107], q[108];
sx q[107];
cx q[108], q[112];
rz(0.33983693) q[107];
sx q[107];
cx q[106], q[107];
sx q[106];
cx q[107], q[108];
rz(0.36136713) q[106];
sx q[106];
cx q[93], q[106];
sx q[93];
cx q[106], q[107];
rz(0.38759673) q[93];
sx q[93];
cx q[87], q[93];
sx q[87];
cx q[93], q[106];
rz(0.42053433) q[87];
sx q[87];
cx q[88], q[87];
sx q[88];
cx q[87], q[93];
rz(0.46364763) q[88];
sx q[88];
cx q[89], q[88];
sx q[89];
cx q[88], q[87];
rz(pi/6) q[89];
sx q[89];
cx q[74], q[89];
sx q[74];
cx q[89], q[88];
rz(0.61547971) q[74];
sx q[74];
cx q[70], q[74];
sx q[70];
cx q[74], q[89];
rz(pi/4) q[70];
sx q[70];
cx q[70], q[74];

// measurement
measure q[70]->c[0];
measure q[74]->c[1];
measure q[89]->c[2];
measure q[88]->c[3];
measure q[87]->c[4];
measure q[93]->c[5];
measure q[106]->c[6];
measure q[107]->c[7];
measure q[108]->c[8];
measure q[112]->c[9];
measure q[126]->c[10];
measure q[125]->c[11];
measure q[124]->c[12];
measure q[123]->c[13];
measure q[122]->c[14];
measure q[121]->c[15];
measure q[120]->c[16];
measure q[119]->c[17];
measure q[118]->c[18];
measure q[110]->c[19];
measure q[100]->c[20];
measure q[101]->c[21];
measure q[102]->c[22];
measure q[92]->c[23];
measure q[83]->c[24];
measure q[84]->c[25];
measure q[85]->c[26];
measure q[73]->c[27];
measure q[66]->c[28];
measure q[65]->c[29];
measure q[64]->c[30];
measure q[63]->c[31];
measure q[62]->c[32];
measure q[61]->c[33];
measure q[60]->c[34];
measure q[53]->c[35];