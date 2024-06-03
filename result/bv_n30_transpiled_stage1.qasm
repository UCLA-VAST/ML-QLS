OPENQASM 2.0;
include "qelib1.inc";
qreg q[127];
creg c[127];
rz(pi/2) q[81];
rz(pi/2) q[7];
rz(pi/2) q[76];
rz(pi/2) q[114];
rz(pi/2) q[79];
rz(pi/2) q[28];
rz(pi/2) q[72];
rz(pi/2) q[21];
rz(pi/2) q[83];
rz(pi/2) q[40];
rz(pi/2) q[92];
rz(pi/2) q[33];
rz(pi/2) q[90];
rz(pi/2) q[85];
rz(pi/2) q[65];
rz(pi/2) q[9];
rz(pi/2) q[116];
rz(pi/2) q[100];
rz(pi/2) q[104];
rz(pi/2) q[80];
sx q[81];
sx q[7];
sx q[76];
sx q[114];
sx q[79];
sx q[28];
sx q[72];
sx q[21];
sx q[83];
sx q[40];
sx q[92];
sx q[33];
sx q[90];
sx q[85];
sx q[65];
sx q[9];
sx q[116];
sx q[100];
sx q[104];
sx q[80];
rz(pi/2) q[81];
rz(pi/2) q[7];
rz(pi/2) q[76];
rz(pi/2) q[114];
rz(pi/2) q[79];
rz(pi/2) q[28];
rz(pi/2) q[72];
rz(pi/2) q[21];
rz(pi/2) q[83];
rz(pi/2) q[40];
rz(pi/2) q[92];
rz(pi/2) q[33];
rz(pi/2) q[90];
rz(pi/2) q[85];
rz(pi/2) q[65];
rz(pi/2) q[9];
rz(pi/2) q[116];
rz(pi/2) q[100];
rz(pi/2) q[104];
rz(-pi/2) q[80];
cx q[81], q[80];
rz(pi/2) q[7];
rz(pi/2) q[40];
rz(pi/2) q[33];
rz(pi/2) q[90];
rz(pi/2) q[65];
rz(pi/2) q[9];
rz(pi/2) q[76];
rz(pi/2) q[116];
rz(pi/2) q[114];
rz(pi/2) q[28];
rz(pi/2) q[21];
cx q[79], q[80];
sx q[7];
sx q[40];
sx q[33];
sx q[90];
sx q[65];
sx q[9];
sx q[76];
sx q[116];
sx q[114];
sx q[28];
sx q[21];
rz(pi/2) q[7];
rz(pi/2) q[40];
rz(pi/2) q[33];
rz(pi/2) q[90];
rz(pi/2) q[65];
rz(pi/2) q[9];
rz(pi/2) q[76];
rz(pi/2) q[116];
rz(pi/2) q[114];
rz(pi/2) q[28];
rz(pi/2) q[21];
swap q[79], q[91];
rz(pi/2) q[79];
rz(pi/2) q[91];
sx q[79];
sx q[91];
rz(pi/2) q[79];
rz(pi/2) q[91];
cx q[79], q[80];
rz(pi/2) q[79];
swap q[80], q[81];
cx q[72], q[81];
rz(pi/2) q[80];
sx q[79];
sx q[80];
rz(pi/2) q[79];
swap q[62], q[72];
rz(pi/2) q[72];
rz(pi/2) q[80];
rz(pi/2) q[62];
sx q[72];
sx q[62];
rz(pi/2) q[72];
rz(pi/2) q[62];
cx q[72], q[81];
rz(pi/2) q[72];
swap q[81], q[82];
rz(pi/2) q[81];
cx q[83], q[82];
sx q[72];
sx q[81];
rz(pi/2) q[72];
rz(pi/2) q[81];
cx q[81], q[82];
rz(pi/2) q[81];
swap q[82], q[83];
cx q[92], q[83];
rz(pi/2) q[82];
sx q[81];
sx q[82];
rz(pi/2) q[81];
swap q[83], q[84];
rz(pi/2) q[83];
rz(pi/2) q[82];
sx q[83];
rz(pi/2) q[83];
cx q[83], q[84];
cx q[85], q[84];
swap q[73], q[85];
rz(pi/2) q[85];
rz(pi/2) q[73];
sx q[85];
sx q[73];
rz(pi/2) q[85];
rz(pi/2) q[73];
cx q[85], q[84];
swap q[85], q[86];
rz(pi/2) q[85];
rz(pi/2) q[86];
sx q[85];
sx q[86];
rz(pi/2) q[85];
rz(pi/2) q[86];
cx q[85], q[84];
rz(pi/2) q[85];
swap q[83], q[84];
rz(pi/2) q[84];
sx q[85];
swap q[83], q[92];
rz(pi/2) q[83];
sx q[84];
rz(pi/2) q[85];
swap q[92], q[102];
rz(pi/2) q[92];
sx q[83];
rz(pi/2) q[84];
swap q[101], q[102];
sx q[92];
cx q[100], q[101];
rz(pi/2) q[83];
rz(pi/2) q[92];
rz(pi/2) q[100];
swap q[101], q[102];
rz(pi/2) q[101];
cx q[92], q[102];
sx q[100];
sx q[101];
rz(pi/2) q[100];
rz(pi/2) q[92];
rz(pi/2) q[101];
sx q[92];
cx q[101], q[102];
rz(pi/2) q[92];
rz(pi/2) q[101];
swap q[102], q[103];
rz(pi/2) q[102];
cx q[104], q[103];
sx q[101];
sx q[102];
rz(pi/2) q[101];
swap q[104], q[105];
rz(pi/2) q[104];
rz(pi/2) q[102];
rz(pi/2) q[105];
sx q[104];
sx q[105];
rz(pi/2) q[104];
rz(pi/2) q[105];
cx q[104], q[103];
cx q[102], q[103];
rz(pi/2) q[104];
sx q[104];
rz(pi/2) q[102];
rz(pi/2) q[104];
sx q[102];
rz(pi/2) q[102];

// measurement
measure q[80]->c[0];
measure q[7]->c[1];
measure q[76]->c[2];
measure q[114]->c[3];
measure q[91]->c[4];
measure q[79]->c[5];
measure q[28]->c[6];
measure q[62]->c[7];
measure q[72]->c[8];
measure q[21]->c[9];
measure q[82]->c[10];
measure q[81]->c[11];
measure q[40]->c[12];
measure q[83]->c[13];
measure q[33]->c[14];
measure q[84]->c[15];
measure q[90]->c[16];
measure q[73]->c[17];
measure q[65]->c[18];
measure q[9]->c[19];
measure q[116]->c[20];
measure q[86]->c[21];
measure q[85]->c[22];
measure q[100]->c[23];
measure q[92]->c[24];
measure q[101]->c[25];
measure q[105]->c[26];
measure q[104]->c[27];
measure q[102]->c[28];
measure q[103]->c[29];
