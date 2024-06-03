OPENQASM 2.0;
include "qelib1.inc";
qreg q[127];
creg c[127];
rz(pi/2) q[77];
rz(-pi) q[75];
rz(-pi) q[92];
rz(-pi) q[85];
rz(-pi) q[80];
rz(-pi) q[71];
rz(-pi) q[62];
rz(-pi) q[93];
rz(-pi) q[65];
rz(-pi) q[58];
rz(-pi) q[42];
rz(-pi) q[61];
rz(-pi) q[53];
rz(-pi) q[79];
rz(-pi) q[98];
rz(-pi) q[101];
rz(-pi) q[76];
rz(-pi) q[84];
rz(-pi) q[83];
rz(-pi) q[81];
rz(-pi) q[78];
rz(-pi) q[72];
rz(-pi) q[87];
rz(-pi) q[66];
rz(-pi) q[59];
rz(-pi) q[41];
rz(-pi) q[60];
rz(-pi) q[63];
rz(-pi) q[82];
rz(-pi) q[99];
rz(-pi) q[100];
sx q[77];
sx q[75];
sx q[92];
sx q[85];
sx q[80];
sx q[71];
sx q[62];
sx q[93];
sx q[65];
sx q[58];
sx q[42];
sx q[61];
sx q[53];
sx q[79];
sx q[98];
sx q[101];
sx q[76];
sx q[84];
sx q[83];
sx q[81];
sx q[78];
sx q[72];
sx q[87];
sx q[66];
sx q[59];
sx q[41];
sx q[60];
sx q[63];
sx q[82];
sx q[99];
sx q[100];
rz(pi/2) q[77];
rz(2.4032558) q[75];
rz(1.0279274) q[92];
rz(2.1498908) q[85];
rz(0.65624515) q[80];
rz(0.16851895) q[71];
rz(3.0258592) q[62];
rz(2.773601) q[93];
rz(0.17562475) q[65];
rz(0.26401255) q[58];
rz(3.0988719) q[42];
rz(1.394253) q[61];
rz(2.3768582) q[53];
rz(0.98853785) q[79];
rz(2.5641275) q[98];
rz(2.8540895) q[101];
rz(0.34614615) q[76];
rz(0.43990665) q[84];
rz(2.4438897) q[83];
rz(1.2955396) q[81];
rz(2.0567011) q[78];
rz(2.5285922) q[72];
rz(1.8206009) q[87];
rz(2.069844) q[66];
rz(1.6367006) q[59];
rz(0.36481265) q[41];
rz(1.1484071) q[60];
rz(0.25950015) q[63];
rz(0.39975555) q[82];
rz(2.1598985) q[99];
rz(0.93872255) q[100];
sx q[75];
sx q[92];
sx q[85];
sx q[80];
sx q[71];
sx q[62];
sx q[93];
sx q[65];
sx q[58];
sx q[42];
sx q[61];
sx q[53];
sx q[79];
sx q[98];
sx q[101];
sx q[76];
sx q[84];
sx q[83];
sx q[81];
sx q[78];
sx q[72];
sx q[87];
sx q[66];
sx q[59];
sx q[41];
sx q[60];
sx q[63];
sx q[82];
sx q[99];
sx q[100];
cx q[75], q[76];
cx q[80], q[81];
cx q[62], q[72];
cx q[93], q[87];
cx q[65], q[66];
cx q[58], q[59];
cx q[42], q[41];
cx q[61], q[60];
cx q[98], q[99];
cx q[101], q[100];
swap q[83], q[84];
rz(pi/2) q[76];
cx q[92], q[83];
cx q[85], q[84];
rz(pi/2) q[81];
rz(pi/2) q[72];
rz(pi/2) q[87];
rz(pi/2) q[66];
rz(pi/2) q[59];
rz(pi/2) q[41];
rz(pi/2) q[60];
rz(pi/2) q[99];
rz(pi/2) q[100];
sx q[76];
rz(pi/2) q[83];
rz(pi/2) q[84];
sx q[81];
sx q[72];
sx q[87];
sx q[66];
sx q[59];
sx q[41];
sx q[60];
sx q[99];
sx q[100];
rz(pi/2) q[76];
sx q[83];
sx q[84];
rz(pi/2) q[81];
rz(pi/2) q[72];
rz(pi/2) q[87];
rz(pi/2) q[66];
rz(pi/2) q[59];
rz(pi/2) q[41];
rz(pi/2) q[60];
rz(pi/2) q[99];
rz(pi/2) q[100];
cx q[75], q[76];
rz(pi/2) q[83];
rz(pi/2) q[84];
cx q[80], q[81];
cx q[62], q[72];
cx q[93], q[87];
cx q[65], q[66];
cx q[58], q[59];
cx q[42], q[41];
cx q[61], q[60];
cx q[98], q[99];
cx q[101], q[100];
rz(-pi/4) q[76];
cx q[92], q[83];
cx q[85], q[84];
rz(-pi/4) q[81];
rz(-pi/4) q[72];
rz(-pi/4) q[87];
rz(-pi/4) q[66];
rz(-pi/4) q[59];
rz(-pi/4) q[41];
rz(-pi/4) q[60];
rz(-pi/4) q[99];
rz(-pi/4) q[100];
swap q[62], q[63];
swap q[79], q[80];
cx q[77], q[76];
rz(-pi/4) q[83];
rz(-pi/4) q[84];
swap q[80], q[81];
swap q[53], q[60];
swap q[98], q[99];
rz(pi/4) q[76];
cx q[81], q[82];
swap q[60], q[61];
cx q[75], q[76];
cx q[61], q[62];
rz(pi/2) q[82];
rz(pi/4) q[75];
rz(-pi/4) q[76];
rz(pi/2) q[62];
sx q[82];
cx q[77], q[76];
sx q[62];
rz(pi/2) q[82];
rz(3*pi/4) q[76];
rz(pi/2) q[62];
cx q[81], q[82];
sx q[76];
cx q[61], q[62];
rz(-pi/4) q[82];
rz(pi/2) q[76];
rz(-pi/4) q[62];
swap q[75], q[76];
cx q[77], q[76];
rz(pi/4) q[77];
rz(-pi/4) q[76];
cx q[77], q[76];
cx q[76], q[75];
swap q[77], q[78];
cx q[71], q[77];
swap q[78], q[79];
rz(pi/2) q[77];
swap q[79], q[80];
sx q[77];
swap q[80], q[81];
rz(pi/2) q[77];
swap q[81], q[82];
cx q[82], q[83];
cx q[71], q[77];
rz(pi/4) q[83];
rz(-pi/4) q[77];
cx q[92], q[83];
rz(-pi/4) q[83];
rz(pi/4) q[92];
cx q[82], q[83];
rz(3*pi/4) q[83];
sx q[83];
rz(pi/2) q[83];
swap q[83], q[92];
cx q[82], q[83];
rz(pi/4) q[82];
rz(-pi/4) q[83];
cx q[82], q[83];
cx q[83], q[92];
swap q[83], q[84];
cx q[82], q[83];
rz(pi/4) q[83];
swap q[83], q[84];
cx q[85], q[84];
rz(-pi/4) q[84];
rz(pi/4) q[85];
swap q[83], q[84];
cx q[82], q[83];
swap q[84], q[85];
rz(3*pi/4) q[83];
sx q[83];
rz(pi/2) q[83];
swap q[83], q[84];
cx q[82], q[83];
rz(pi/4) q[82];
rz(-pi/4) q[83];
cx q[82], q[83];
cx q[83], q[84];
swap q[81], q[82];
swap q[80], q[81];
cx q[80], q[79];
rz(pi/4) q[79];
cx q[78], q[79];
rz(-pi/4) q[79];
rz(pi/4) q[78];
cx q[80], q[79];
rz(3*pi/4) q[79];
sx q[79];
rz(pi/2) q[79];
swap q[79], q[80];
cx q[79], q[78];
rz(pi/4) q[79];
rz(-pi/4) q[78];
cx q[79], q[78];
swap q[77], q[78];
cx q[79], q[78];
swap q[71], q[77];
rz(pi/4) q[78];
cx q[77], q[78];
rz(-pi/4) q[78];
rz(pi/4) q[77];
cx q[79], q[78];
rz(3*pi/4) q[78];
sx q[78];
rz(pi/2) q[78];
swap q[77], q[78];
cx q[79], q[78];
rz(pi/4) q[79];
rz(-pi/4) q[78];
cx q[79], q[78];
cx q[78], q[77];
swap q[79], q[80];
swap q[80], q[81];
swap q[71], q[77];
swap q[78], q[79];
cx q[77], q[78];
cx q[81], q[72];
rz(pi/4) q[72];
swap q[62], q[72];
cx q[63], q[62];
swap q[72], q[81];
rz(-pi/4) q[62];
rz(pi/4) q[63];
cx q[72], q[62];
rz(3*pi/4) q[62];
sx q[62];
rz(pi/2) q[62];
swap q[62], q[63];
cx q[72], q[62];
rz(pi/4) q[72];
rz(-pi/4) q[62];
cx q[72], q[62];
cx q[62], q[63];
swap q[72], q[81];
swap q[81], q[82];
swap q[82], q[83];
swap q[83], q[84];
swap q[84], q[85];
swap q[85], q[86];
cx q[86], q[87];
rz(pi/4) q[87];
cx q[93], q[87];
rz(-pi/4) q[87];
rz(pi/4) q[93];
cx q[86], q[87];
rz(3*pi/4) q[87];
sx q[87];
rz(pi/2) q[87];
swap q[87], q[93];
cx q[86], q[87];
rz(pi/4) q[86];
rz(-pi/4) q[87];
cx q[86], q[87];
cx q[87], q[93];
swap q[85], q[86];
swap q[73], q[85];
cx q[73], q[66];
rz(pi/4) q[66];
cx q[65], q[66];
rz(-pi/4) q[66];
rz(pi/4) q[65];
cx q[73], q[66];
rz(3*pi/4) q[66];
sx q[66];
rz(pi/2) q[66];
swap q[66], q[73];
cx q[66], q[65];
rz(pi/4) q[66];
rz(-pi/4) q[65];
cx q[66], q[65];
swap q[65], q[66];
cx q[66], q[73];
swap q[64], q[65];
swap q[63], q[64];
swap q[62], q[63];
swap q[61], q[62];
swap q[60], q[61];
swap q[62], q[72];
cx q[60], q[59];
rz(pi/4) q[59];
cx q[58], q[59];
rz(-pi/4) q[59];
rz(pi/4) q[58];
cx q[60], q[59];
rz(3*pi/4) q[59];
sx q[59];
rz(pi/2) q[59];
swap q[58], q[59];
cx q[60], q[59];
rz(pi/4) q[60];
rz(-pi/4) q[59];
cx q[60], q[59];
cx q[59], q[58];
swap q[53], q[60];
cx q[53], q[41];
rz(pi/4) q[41];
cx q[42], q[41];
rz(pi/4) q[42];
rz(-pi/4) q[41];
cx q[53], q[41];
rz(3*pi/4) q[41];
sx q[41];
rz(pi/2) q[41];
swap q[41], q[42];
cx q[53], q[41];
rz(pi/4) q[53];
rz(-pi/4) q[41];
cx q[53], q[41];
cx q[41], q[42];
cx q[53], q[60];
rz(pi/4) q[60];
cx q[61], q[60];
rz(pi/4) q[61];
rz(-pi/4) q[60];
cx q[53], q[60];
rz(3*pi/4) q[60];
sx q[60];
rz(pi/2) q[60];
swap q[53], q[60];
cx q[60], q[61];
rz(pi/4) q[60];
rz(-pi/4) q[61];
cx q[60], q[61];
swap q[60], q[61];
cx q[60], q[53];
cx q[61], q[62];
rz(pi/4) q[62];
cx q[72], q[62];
rz(pi/4) q[72];
rz(-pi/4) q[62];
cx q[61], q[62];
rz(3*pi/4) q[62];
sx q[62];
rz(pi/2) q[62];
swap q[61], q[62];
cx q[62], q[72];
rz(pi/4) q[62];
rz(-pi/4) q[72];
cx q[62], q[72];
swap q[62], q[72];
cx q[62], q[61];
cx q[72], q[81];
rz(pi/4) q[81];
cx q[80], q[81];
rz(pi/4) q[80];
rz(-pi/4) q[81];
cx q[72], q[81];
rz(3*pi/4) q[81];
sx q[81];
rz(pi/2) q[81];
swap q[72], q[81];
cx q[81], q[80];
rz(pi/4) q[81];
rz(-pi/4) q[80];
cx q[81], q[80];
swap q[80], q[81];
cx q[81], q[72];
swap q[79], q[80];
swap q[79], q[91];
cx q[91], q[98];
rz(pi/4) q[98];
cx q[99], q[98];
rz(pi/4) q[99];
rz(-pi/4) q[98];
cx q[91], q[98];
rz(3*pi/4) q[98];
sx q[98];
rz(pi/2) q[98];
swap q[91], q[98];
cx q[98], q[99];
rz(pi/4) q[98];
rz(-pi/4) q[99];
cx q[98], q[99];
swap q[98], q[99];
cx q[98], q[91];
cx q[99], q[100];
rz(pi/4) q[100];
cx q[101], q[100];
rz(pi/4) q[101];
rz(-pi/4) q[100];
cx q[99], q[100];
rz(3*pi/4) q[100];
sx q[100];
rz(pi/2) q[100];
swap q[100], q[101];
cx q[99], q[100];
rz(pi/4) q[99];
rz(-pi/4) q[100];
cx q[99], q[100];
rz(pi/2) q[99];
cx q[100], q[101];
sx q[99];
rz(pi/2) q[99];

// measurement
measure q[99]->c[0];
measure q[76]->c[1];
measure q[84]->c[2];
measure q[82]->c[3];
measure q[77]->c[4];
measure q[80]->c[5];
measure q[63]->c[6];
measure q[87]->c[7];
measure q[66]->c[8];
measure q[59]->c[9];
measure q[41]->c[10];
measure q[60]->c[11];
measure q[62]->c[12];
measure q[81]->c[13];
measure q[98]->c[14];
measure q[100]->c[15];
measure q[75]->c[16];
measure q[92]->c[17];
measure q[83]->c[18];
measure q[78]->c[19];
measure q[71]->c[20];
measure q[64]->c[21];
measure q[93]->c[22];
measure q[73]->c[23];
measure q[58]->c[24];
measure q[42]->c[25];
measure q[53]->c[26];
measure q[61]->c[27];
measure q[72]->c[28];
measure q[91]->c[29];
measure q[101]->c[30];
