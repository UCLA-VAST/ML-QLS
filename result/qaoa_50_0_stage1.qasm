OPENQASM 2.0;
include "qelib1.inc";
qreg q[127];
creg c[127];
rzz(pi/4) q[51], q[36];
rzz(pi/4) q[81], q[82];
rzz(pi/4) q[52], q[37];
rzz(pi/4) q[44], q[43];
rzz(pi/4) q[66], q[67];
rzz(pi/4) q[21], q[20];
rzz(pi/4) q[46], q[45];
rzz(pi/4) q[49], q[48];
rzz(pi/4) q[27], q[26];
rzz(pi/4) q[60], q[61];
swap q[25], q[26];
swap q[27], q[28];
swap q[46], q[47];
swap q[34], q[43];
swap q[44], q[45];
swap q[52], q[56];
swap q[37], q[38];
swap q[67], q[68];
swap q[65], q[66];
swap q[59], q[60];
swap q[61], q[62];
swap q[32], q[36];
swap q[50], q[51];
swap q[72], q[81];
rzz(pi/4) q[81], q[72];
rzz(pi/4) q[26], q[27];
rzz(pi/4) q[53], q[60];
swap q[24], q[25];
swap q[28], q[29];
swap q[43], q[44];
swap q[56], q[57];
swap q[45], q[54];
swap q[38], q[39];
swap q[55], q[68];
swap q[31], q[32];
swap q[47], q[48];
rzz(pi/4) q[29], q[30];
rzz(pi/4) q[57], q[58];
rzz(pi/4) q[68], q[67];
rzz(pi/4) q[81], q[82];
rzz(pi/4) q[26], q[25];
swap q[23], q[24];
swap q[42], q[43];
swap q[33], q[39];
swap q[54], q[64];
swap q[41], q[53];
swap q[27], q[28];
swap q[49], q[55];
swap q[22], q[23];
swap q[58], q[71];
swap q[24], q[34];
swap q[20], q[33];
swap q[64], q[65];
swap q[41], q[42];
swap q[16], q[26];
swap q[49], q[50];
swap q[30], q[31];
swap q[66], q[67];
swap q[55], q[68];
rzz(pi/4) q[32], q[31];
rzz(pi/4) q[41], q[53];
rzz(pi/4) q[68], q[67];
rzz(pi/4) q[43], q[34];
rzz(pi/4) q[59], q[58];
swap q[71], q[77];
swap q[23], q[24];
swap q[33], q[39];
swap q[63], q[64];
swap q[21], q[22];
rzz(pi/4) q[64], q[54];
swap q[77], q[78];
swap q[24], q[34];
swap q[42], q[43];
swap q[62], q[63];
swap q[22], q[23];
swap q[20], q[21];
swap q[59], q[60];
swap q[40], q[41];
swap q[32], q[36];
rzz(pi/4) q[36], q[51];
rzz(pi/4) q[20], q[33];
rzz(pi/4) q[22], q[21];
rzz(pi/4) q[39], q[40];
rzz(pi/4) q[63], q[64];
rzz(pi/4) q[59], q[58];
swap q[78], q[79];
swap q[34], q[43];
swap q[45], q[54];
swap q[24], q[25];
swap q[53], q[60];
swap q[62], q[72];
rzz(pi/4) q[53], q[41];
rzz(pi/4) q[61], q[62];
rzz(pi/4) q[58], q[57];
rzz(pi/4) q[21], q[20];
rzz(pi/4) q[79], q[80];
rzz(pi/4) q[25], q[26];
rzz(pi/4) q[59], q[60];
rzz(pi/4) q[23], q[24];
swap q[45], q[46];
swap q[54], q[64];
swap q[50], q[51];
swap q[33], q[39];
rzz(pi/4) q[41], q[42];
rzz(pi/4) q[33], q[39];
rzz(pi/4) q[61], q[60];
rzz(pi/4) q[64], q[63];
swap q[46], q[47];
swap q[80], q[81];
swap q[22], q[23];
rzz(pi/4) q[46], q[45];
rzz(pi/4) q[65], q[64];
rzz(pi/4) q[72], q[81];
rzz(pi/4) q[22], q[23];
rzz(pi/4) q[80], q[79];
swap q[35], q[47];
swap q[40], q[41];
rzz(pi/4) q[40], q[39];
rzz(pi/4) q[48], q[47];
rzz(pi/4) q[81], q[82];
swap q[28], q[35];
swap q[44], q[45];
swap q[62], q[72];
rzz(pi/4) q[43], q[44];
rzz(pi/4) q[62], q[61];
rzz(pi/4) q[49], q[48];
rzz(pi/4) q[27], q[28];
swap q[46], q[47];
rzz(pi/4) q[45], q[46];
swap q[28], q[35];
swap q[26], q[27];
swap q[34], q[43];
swap q[49], q[50];
rzz(pi/4) q[45], q[44];
rzz(pi/4) q[50], q[51];
rzz(pi/4) q[47], q[35];
rzz(pi/4) q[26], q[16];
rzz(pi/4) q[34], q[24];
rzz(pi/4) q[49], q[55];
swap q[28], q[29];
rzz(pi/4) q[36], q[51];
rzz(pi/4) q[26], q[25];
rzz(pi/4) q[43], q[34];
swap q[27], q[28];
swap q[55], q[68];
swap q[45], q[54];
rzz(pi/4) q[64], q[54];
rzz(pi/4) q[45], q[46];
rzz(pi/4) q[55], q[49];
rzz(pi/4) q[27], q[28];
rzz(pi/4) q[67], q[68];
rzz(pi/4) q[43], q[42];
rzz(pi/4) q[67], q[66];
rzz(pi/4) q[29], q[28];
rzz(pi/4) q[30], q[29];
rzz(pi/4) q[65], q[66];
rzz(pi/4) q[30], q[31];

// measurement
measure q[36]->c[0];
measure q[65]->c[1];
measure q[64]->c[2];
measure q[55]->c[3];
measure q[43]->c[4];
measure q[34]->c[5];
measure q[53]->c[6];
measure q[80]->c[7];
measure q[59]->c[8];
measure q[63]->c[9];
measure q[58]->c[10];
measure q[40]->c[11];
measure q[22]->c[12];
measure q[26]->c[13];
measure q[16]->c[14];
measure q[33]->c[15];
measure q[57]->c[16];
measure q[23]->c[17];
measure q[50]->c[18];
measure q[47]->c[19];
measure q[42]->c[20];
measure q[25]->c[21];
measure q[62]->c[22];
measure q[67]->c[23];
measure q[30]->c[24];
measure q[54]->c[25];
measure q[21]->c[26];
measure q[49]->c[27];
measure q[27]->c[28];
measure q[44]->c[29];
measure q[31]->c[30];
measure q[20]->c[31];
measure q[45]->c[32];
measure q[79]->c[33];
measure q[48]->c[34];
measure q[51]->c[35];
measure q[41]->c[36];
measure q[29]->c[37];
measure q[35]->c[38];
measure q[68]->c[39];
measure q[66]->c[40];
measure q[61]->c[41];
measure q[39]->c[42];
measure q[72]->c[43];
measure q[46]->c[44];
measure q[24]->c[45];
measure q[81]->c[46];
measure q[60]->c[47];
measure q[82]->c[48];
measure q[28]->c[49];
