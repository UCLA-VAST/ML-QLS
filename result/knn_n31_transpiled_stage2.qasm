OPENQASM 2.0;
include "qelib1.inc";
qreg q[36];
creg c[36];
rz(pi/2) q[11];
rz(-pi) q[9];
rz(-pi) q[3];
rz(-pi) q[1];
rz(-pi) q[19];
rz(-pi) q[18];
rz(-pi) q[0];
rz(-pi) q[31];
rz(-pi) q[32];
rz(-pi) q[22];
rz(-pi) q[17];
rz(-pi) q[10];
rz(-pi) q[4];
rz(-pi) q[2];
rz(-pi) q[20];
rz(-pi) q[7];
rz(-pi) q[12];
rz(-pi) q[6];
rz(-pi) q[30];
rz(-pi) q[33];
rz(-pi) q[23];
rz(-pi) q[35];
rz(-pi) q[15];
sx q[11];
sx q[9];
sx q[3];
sx q[1];
sx q[19];
sx q[18];
sx q[0];
sx q[31];
sx q[32];
sx q[22];
sx q[17];
sx q[10];
sx q[4];
sx q[2];
sx q[20];
sx q[7];
sx q[12];
sx q[6];
sx q[30];
sx q[33];
sx q[23];
sx q[35];
sx q[15];
rz(pi/2) q[11];
rz(2.4032558) q[9];
rz(1.0279274) q[3];
rz(2.1498908) q[1];
rz(0.65624515) q[19];
rz(3.0258592) q[18];
rz(2.773601) q[0];
rz(0.17562475) q[31];
rz(3.0988719) q[32];
rz(2.3768582) q[22];
rz(0.98853785) q[17];
rz(0.34614615) q[10];
rz(0.43990665) q[4];
rz(2.4438897) q[2];
rz(1.2955396) q[20];
rz(2.0567011) q[7];
rz(2.5285922) q[12];
rz(1.8206009) q[6];
rz(2.069844) q[30];
rz(0.36481265) q[33];
rz(0.39975555) q[23];
rz(2.1598985) q[35];
rz(0.93872255) q[15];
sx q[9];
sx q[3];
sx q[1];
sx q[19];
sx q[18];
sx q[0];
sx q[31];
sx q[32];
sx q[22];
sx q[17];
sx q[10];
sx q[4];
sx q[2];
sx q[20];
sx q[7];
sx q[12];
sx q[6];
sx q[30];
sx q[33];
sx q[23];
sx q[35];
sx q[15];
cx q[9], q[10];
cx q[3], q[4];
cx q[1], q[2];
cx q[19], q[20];
cx q[18], q[12];
cx q[0], q[6];
cx q[31], q[30];
cx q[32], q[33];
cx q[17], q[23];
rz(pi/2) q[10];
rz(pi/2) q[4];
rz(pi/2) q[2];
rz(pi/2) q[20];
rz(pi/2) q[12];
rz(pi/2) q[6];
rz(pi/2) q[30];
rz(pi/2) q[33];
rz(pi/2) q[23];
sx q[10];
sx q[4];
sx q[2];
sx q[20];
sx q[12];
sx q[6];
sx q[30];
sx q[33];
sx q[23];
rz(pi/2) q[10];
rz(pi/2) q[4];
rz(pi/2) q[2];
rz(pi/2) q[20];
rz(pi/2) q[12];
rz(pi/2) q[6];
rz(pi/2) q[30];
rz(pi/2) q[33];
rz(pi/2) q[23];
cx q[9], q[10];
cx q[3], q[4];
cx q[1], q[2];
cx q[19], q[20];
cx q[18], q[12];
cx q[0], q[6];
cx q[31], q[30];
cx q[32], q[33];
cx q[17], q[23];
rz(-pi/4) q[10];
rz(-pi/4) q[4];
rz(-pi/4) q[2];
rz(-pi/4) q[20];
rz(-pi/4) q[12];
rz(-pi/4) q[6];
rz(-pi/4) q[30];
rz(-pi/4) q[33];
rz(-pi/4) q[23];
cx q[11], q[10];
rz(pi/4) q[10];
cx q[9], q[10];
rz(pi/4) q[9];
swap q[10], q[11];
rz(-pi/4) q[11];
cx q[10], q[11];
cx q[10], q[9];
rz(3*pi/4) q[11];
rz(pi/4) q[10];
rz(-pi/4) q[9];
sx q[11];
cx q[10], q[9];
rz(pi/2) q[11];
cx q[10], q[4];
rz(pi/4) q[4];
cx q[3], q[4];
rz(-pi/4) q[4];
rz(pi/4) q[3];
cx q[10], q[4];
rz(3*pi/4) q[4];
swap q[9], q[10];
cx q[10], q[11];
sx q[4];
cx q[9], q[3];
rz(pi/2) q[4];
rz(pi/4) q[9];
rz(-pi/4) q[3];
cx q[9], q[3];
cx q[3], q[4];
swap q[8], q[9];
cx q[8], q[2];
rz(pi/4) q[2];
cx q[1], q[2];
rz(-pi/4) q[2];
cx q[8], q[2];
swap q[1], q[2];
rz(3*pi/4) q[1];
rz(pi/4) q[2];
sx q[1];
cx q[8], q[2];
rz(pi/2) q[1];
rz(pi/4) q[8];
rz(-pi/4) q[2];
cx q[8], q[2];
cx q[2], q[1];
swap q[8], q[14];
rz(-pi) q[8];
cx q[14], q[20];
sx q[8];
rz(pi/4) q[20];
rz(0.16851895) q[8];
cx q[19], q[20];
sx q[8];
rz(-pi/4) q[20];
rz(pi/4) q[19];
cx q[14], q[20];
cx q[8], q[7];
rz(3*pi/4) q[20];
rz(pi/2) q[7];
swap q[13], q[14];
rz(-pi) q[14];
sx q[20];
cx q[13], q[19];
sx q[7];
sx q[14];
rz(pi/2) q[20];
rz(pi/4) q[13];
rz(-pi/4) q[19];
rz(pi/2) q[7];
rz(2.8540895) q[14];
cx q[13], q[19];
cx q[8], q[7];
sx q[14];
cx q[19], q[20];
rz(-pi/4) q[7];
cx q[13], q[7];
cx q[14], q[15];
rz(pi/4) q[7];
rz(pi/2) q[15];
cx q[8], q[7];
sx q[15];
rz(-pi/4) q[7];
rz(pi/2) q[15];
cx q[13], q[7];
cx q[14], q[15];
rz(-pi/4) q[15];
swap q[7], q[8];
rz(3*pi/4) q[8];
rz(pi/4) q[7];
sx q[8];
cx q[13], q[7];
rz(pi/2) q[8];
rz(pi/4) q[13];
rz(-pi/4) q[7];
cx q[13], q[7];
cx q[7], q[8];
cx q[13], q[12];
rz(pi/4) q[12];
cx q[18], q[12];
rz(pi/4) q[18];
swap q[12], q[13];
rz(-pi/4) q[13];
cx q[12], q[13];
rz(3*pi/4) q[13];
cx q[12], q[18];
sx q[13];
rz(pi/4) q[12];
rz(-pi/4) q[18];
rz(pi/2) q[13];
cx q[12], q[18];
cx q[12], q[6];
rz(pi/4) q[6];
cx q[0], q[6];
rz(-pi/4) q[6];
cx q[12], q[6];
swap q[0], q[6];
rz(3*pi/4) q[0];
rz(pi/4) q[6];
sx q[0];
cx q[12], q[6];
rz(pi/2) q[0];
rz(pi/4) q[12];
rz(-pi/4) q[6];
cx q[12], q[6];
cx q[6], q[0];
swap q[12], q[18];
cx q[12], q[13];
swap q[18], q[24];
cx q[24], q[30];
rz(pi/4) q[30];
cx q[31], q[30];
rz(-pi/4) q[30];
rz(pi/4) q[31];
cx q[24], q[30];
rz(3*pi/4) q[30];
swap q[24], q[25];
rz(-pi) q[24];
sx q[30];
cx q[25], q[31];
sx q[24];
rz(pi/2) q[30];
rz(pi/4) q[25];
rz(-pi/4) q[31];
rz(0.26401255) q[24];
cx q[25], q[31];
sx q[24];
cx q[31], q[30];
swap q[25], q[26];
rz(-pi) q[25];
sx q[25];
rz(1.6367006) q[25];
sx q[25];
cx q[24], q[25];
rz(pi/2) q[25];
sx q[25];
rz(pi/2) q[25];
cx q[24], q[25];
rz(-pi/4) q[25];
cx q[26], q[25];
rz(pi/4) q[25];
cx q[24], q[25];
rz(-pi/4) q[25];
cx q[26], q[25];
swap q[24], q[25];
rz(3*pi/4) q[24];
rz(pi/4) q[25];
sx q[24];
cx q[26], q[25];
rz(pi/2) q[24];
rz(pi/4) q[26];
rz(-pi/4) q[25];
cx q[26], q[25];
cx q[25], q[24];
swap q[26], q[27];
cx q[27], q[33];
rz(pi/4) q[33];
cx q[32], q[33];
rz(pi/4) q[32];
swap q[27], q[33];
rz(-pi/4) q[27];
cx q[33], q[27];
cx q[33], q[32];
rz(3*pi/4) q[27];
rz(pi/4) q[33];
rz(-pi/4) q[32];
sx q[27];
cx q[33], q[32];
rz(pi/2) q[27];
swap q[33], q[34];
swap q[26], q[32];
rz(-pi) q[32];
cx q[26], q[27];
rz(-pi) q[33];
sx q[32];
sx q[33];
swap q[21], q[27];
rz(1.394253) q[32];
rz(1.1484071) q[33];
sx q[32];
sx q[33];
cx q[32], q[33];
rz(pi/2) q[33];
sx q[33];
rz(pi/2) q[33];
cx q[32], q[33];
rz(-pi/4) q[33];
cx q[34], q[33];
rz(pi/4) q[33];
cx q[32], q[33];
rz(-pi/4) q[33];
cx q[34], q[33];
swap q[32], q[33];
rz(pi/4) q[33];
rz(3*pi/4) q[32];
cx q[34], q[33];
sx q[32];
rz(pi/4) q[34];
rz(-pi/4) q[33];
rz(pi/2) q[32];
cx q[34], q[33];
cx q[33], q[32];
swap q[28], q[34];
rz(-pi) q[34];
swap q[28], q[29];
sx q[34];
rz(-pi) q[28];
rz(2.5641275) q[34];
sx q[28];
sx q[34];
rz(0.25950015) q[28];
sx q[28];
cx q[34], q[35];
cx q[22], q[28];
rz(pi/2) q[35];
rz(pi/2) q[28];
sx q[35];
sx q[28];
rz(pi/2) q[35];
rz(pi/2) q[28];
cx q[34], q[35];
cx q[22], q[28];
rz(-pi/4) q[35];
rz(-pi/4) q[28];
cx q[29], q[28];
rz(pi/4) q[28];
cx q[22], q[28];
rz(-pi/4) q[28];
cx q[29], q[28];
swap q[22], q[28];
rz(pi/4) q[28];
cx q[29], q[28];
rz(pi/4) q[29];
rz(-pi/4) q[28];
cx q[29], q[28];
cx q[29], q[23];
rz(pi/4) q[23];
cx q[17], q[23];
rz(-pi/4) q[23];
cx q[29], q[23];
swap q[17], q[23];
rz(pi/4) q[23];
rz(3*pi/4) q[17];
cx q[29], q[23];
sx q[17];
rz(pi/4) q[29];
rz(-pi/4) q[23];
rz(pi/2) q[17];
cx q[29], q[23];
cx q[23], q[17];
cx q[29], q[35];
rz(pi/4) q[35];
cx q[34], q[35];
rz(pi/4) q[34];
rz(-pi/4) q[35];
cx q[29], q[35];
rz(3*pi/4) q[35];
swap q[28], q[29];
cx q[28], q[34];
sx q[35];
rz(pi/4) q[28];
rz(-pi/4) q[34];
rz(pi/2) q[35];
cx q[28], q[34];
cx q[34], q[35];
swap q[22], q[28];
rz(3*pi/4) q[28];
swap q[21], q[22];
sx q[28];
cx q[21], q[15];
rz(pi/2) q[28];
rz(pi/4) q[15];
cx q[29], q[28];
cx q[14], q[15];
rz(-pi/4) q[15];
cx q[21], q[15];
swap q[14], q[15];
rz(pi/4) q[15];
rz(3*pi/4) q[14];
cx q[21], q[15];
sx q[14];
rz(pi/4) q[21];
rz(-pi/4) q[15];
rz(pi/2) q[14];
cx q[21], q[15];
rz(pi/2) q[21];
cx q[15], q[14];
sx q[21];
rz(pi/2) q[21];

// measurement
measure q[21]->c[0];
measure q[10]->c[1];
measure q[3]->c[2];
measure q[2]->c[3];
measure q[19]->c[4];
measure q[7]->c[5];
measure q[12]->c[6];
measure q[6]->c[7];
measure q[31]->c[8];
measure q[25]->c[9];
measure q[26]->c[10];
measure q[33]->c[11];
measure q[29]->c[12];
measure q[23]->c[13];
measure q[34]->c[14];
measure q[15]->c[15];
measure q[11]->c[16];
measure q[4]->c[17];
measure q[1]->c[18];
measure q[20]->c[19];
measure q[8]->c[20];
measure q[13]->c[21];
measure q[0]->c[22];
measure q[30]->c[23];
measure q[24]->c[24];
measure q[22]->c[25];
measure q[32]->c[26];
measure q[28]->c[27];
measure q[17]->c[28];
measure q[35]->c[29];
measure q[14]->c[30];