OPENQASM 2.0;
include "qelib1.inc";
qreg q[54];
creg c[54];
x q[43];
x q[12];
x q[47];
x q[36];
x q[0];
x q[21];
x q[2];
x q[20];
x q[11];
x q[7];
x q[29];
x q[49];
x q[24];
x q[28];
x q[23];
x q[48];
x q[50];
x q[15];
x q[22];
x q[4];
x q[30];
x q[34];
x q[17];
x q[18];
x q[52];
x q[14];
x q[13];
x q[8];
x q[44];
x q[3];
cx q[33], q[39];
cx q[9], q[16];
cx q[35], q[41];
cx q[45], q[51];
cx q[37], q[42];
cx q[10], q[5];
cx q[53], q[46];
cx q[6], q[1];
cx q[26], q[19];
cx q[27], q[32];
cx q[31], q[38];
x q[40];
x q[34];
x q[44];
x q[49];
x q[33];
x q[13];
x q[21];
x q[24];
x q[47];
x q[42];
x q[45];
x q[6];
x q[10];
x q[8];
x q[18];
x q[35];
x q[31];
x q[16];
x q[28];
x q[5];
x q[51];
x q[27];
x q[3];
x q[0];
x q[50];
x q[23];
x q[22];
x q[36];
x q[4];
x q[29];
cx q[37], q[43];
cx q[25], q[30];
cx q[46], q[52];
cx q[20], q[26];
cx q[39], q[32];
cx q[14], q[19];
cx q[7], q[1];
cx q[15], q[9];
x q[40];
x q[2];
x q[37];
x q[51];
x q[45];
x q[6];
x q[16];
x q[13];
x q[23];
x q[34];
x q[0];
x q[36];
x q[30];
x q[26];
x q[49];
x q[21];
x q[29];
x q[39];
x q[1];
x q[35];
x q[25];
x q[8];
x q[52];
x q[24];
x q[33];
x q[28];
x q[44];
cx q[12], q[18];
cx q[7], q[14];
cx q[31], q[38];
cx q[53], q[47];
cx q[42], q[48];
cx q[43], q[50];
cx q[22], q[17];
cx q[10], q[4];
cx q[9], q[3];
cx q[15], q[20];
cx q[27], q[32];
cx q[11], q[5];
x q[19];
cx q[46], q[41];
x q[37];
x q[53];
x q[47];
x q[12];
x q[3];
x q[49];
x q[0];
x q[20];
x q[43];
x q[27];
x q[24];
x q[30];
x q[33];
x q[42];
x q[5];
x q[11];
x q[48];
x q[52];
x q[14];
x q[36];
x q[9];
cx q[22], q[28];
cx q[39], q[44];
cx q[45], q[51];
cx q[40], q[34];
cx q[25], q[31];
cx q[32], q[26];
cx q[35], q[29];
cx q[10], q[4];
cx q[2], q[8];
cx q[17], q[23];
cx q[6], q[1];
cx q[15], q[21];
cx q[7], q[13];
x q[18];
x q[46];
x q[19];
x q[50];
x q[16];
x q[37];
x q[34];
x q[40];
x q[52];
x q[43];
x q[1];
x q[7];
x q[0];
x q[27];
x q[53];
x q[39];
x q[30];
x q[2];
x q[23];
x q[12];
x q[9];
x q[29];
x q[22];
x q[14];
x q[49];
x q[3];
x q[45];
x q[8];
x q[47];
cx q[44], q[38];
cx q[6], q[13];
cx q[15], q[21];
cx q[25], q[31];
cx q[35], q[41];
cx q[42], q[48];
cx q[33], q[28];
cx q[11], q[5];
cx q[32], q[26];
x q[51];
x q[10];
x q[19];
x q[20];
x q[46];
cx q[24], q[18];
x q[16];
x q[37];
x q[34];
x q[43];
x q[1];
x q[35];
x q[8];
x q[25];
x q[21];
x q[49];
x q[15];
x q[0];
x q[7];
x q[27];
x q[48];
x q[41];
x q[2];
x q[26];
x q[32];
x q[12];
x q[14];
x q[42];
cx q[6], q[13];
cx q[50], q[44];
cx q[31], q[38];
cx q[33], q[39];
cx q[11], q[5];
cx q[9], q[4];
cx q[30], q[36];
cx q[45], q[52];
cx q[53], q[47];
cx q[22], q[17];
cx q[29], q[23];
x q[18];
x q[24];
x q[28];
x q[46];
x q[51];
x q[16];
x q[31];
x q[22];
x q[2];
x q[34];
x q[27];
x q[0];
x q[41];
x q[17];
x q[44];
x q[35];
x q[45];
x q[7];
x q[50];
x q[48];
x q[38];
cx q[37], q[30];
cx q[39], q[32];
cx q[29], q[23];
cx q[33], q[40];
cx q[36], q[42];
cx q[14], q[19];
cx q[15], q[9];
cx q[49], q[43];
cx q[20], q[26];
cx q[11], q[5];
cx q[53], q[47];
cx q[6], q[13];
cx q[8], q[3];
cx q[10], q[4];
x q[18];
x q[12];
x q[24];
x q[25];
x q[51];
x q[46];
x q[21];
x q[52];
x q[40];
x q[30];
x q[38];
x q[2];
x q[45];
x q[27];
x q[3];
x q[11];
x q[5];
x q[4];
x q[19];
x q[15];
x q[9];
x q[48];
x q[17];
x q[8];
cx q[37], q[31];
cx q[7], q[1];
cx q[36], q[42];
cx q[6], q[0];
cx q[34], q[41];
cx q[32], q[26];
cx q[35], q[29];
cx q[22], q[28];
cx q[20], q[14];
cx q[49], q[43];
cx q[50], q[44];
cx q[33], q[39];
cx q[53], q[47];
x q[24];
x q[25];
x q[46];
x q[51];
x q[13];
x q[10];
cx q[16], q[21];
x q[18];
x q[52];
x q[31];
x q[3];
x q[20];
x q[14];
x q[22];
x q[33];
x q[1];
x q[38];
x q[4];
x q[0];
x q[35];
x q[9];
x q[17];
x q[48];
x q[53];
x q[2];
x q[50];
cx q[26], q[19];
cx q[30], q[36];
cx q[49], q[43];
cx q[6], q[12];
cx q[39], q[32];
cx q[37], q[42];
cx q[47], q[41];
cx q[34], q[28];
cx q[29], q[23];
cx q[11], q[5];
cx q[40], q[45];
cx q[15], q[8];
x q[51];
x q[7];
x q[13];
x q[24];
cx q[16], q[21];
x q[18];
x q[52];
x q[31];
x q[26];
x q[40];
x q[20];
x q[3];
x q[32];
x q[2];
x q[8];
x q[23];
x q[35];
x q[30];
x q[41];
x q[11];
x q[22];
x q[39];
x q[15];
x q[45];
x q[12];
x q[38];
x q[14];
x q[47];
x q[29];
cx q[34], q[28];
cx q[37], q[43];
cx q[10], q[5];
cx q[9], q[4];
cx q[50], q[44];
cx q[25], q[19];
cx q[49], q[42];
cx q[33], q[27];
cx q[6], q[1];
x q[17];
x q[51];
x q[24];
x q[21];
x q[36];
x q[48];
cx q[7], q[13];
cx q[53], q[46];
x q[18];
x q[52];
x q[31];
x q[33];
x q[35];
x q[5];
x q[2];
x q[22];
x q[11];
x q[1];
x q[12];
x q[43];
x q[47];
x q[15];
x q[37];
x q[20];
x q[23];
x q[29];
x q[41];
cx q[25], q[30];
cx q[10], q[16];
cx q[6], q[0];
cx q[50], q[44];
cx q[26], q[19];
cx q[14], q[8];
cx q[34], q[28];
cx q[49], q[42];
cx q[40], q[45];
cx q[39], q[32];
cx q[9], q[4];
x q[53];
x q[17];
x q[51];
x q[36];
x q[24];
x q[46];
x q[52];
x q[23];
x q[39];
x q[42];
x q[45];
x q[25];
x q[0];
x q[12];
x q[40];
x q[30];
x q[15];
x q[35];
x q[49];
x q[32];
x q[29];
x q[50];
x q[20];
x q[11];
x q[5];
x q[47];
cx q[31], q[38];
cx q[34], q[41];
cx q[7], q[14];
cx q[9], q[3];
cx q[16], q[21];
cx q[37], q[43];
cx q[2], q[8];
cx q[6], q[1];
cx q[10], q[4];
cx q[33], q[27];
cx q[13], q[19];
cx q[22], q[28];
x q[36];
x q[44];
x q[51];
x q[24];
x q[53];
x q[52];
x q[38];
x q[3];
x q[45];
x q[28];
x q[40];
x q[1];
x q[6];
x q[23];
x q[37];
x q[15];
x q[5];
x q[0];
x q[33];
x q[39];
x q[31];
x q[2];
x q[49];
x q[11];
x q[29];
x q[47];
x q[20];
x q[21];
cx q[14], q[19];
cx q[7], q[13];
cx q[12], q[18];
cx q[10], q[17];
cx q[42], q[48];
cx q[25], q[30];
cx q[35], q[41];
cx q[43], q[50];
cx q[32], q[26];
cx q[9], q[4];
x q[8];
x q[44];
x q[24];
x q[22];
x q[16];
x q[51];
x q[52];
x q[38];
x q[10];
x q[33];
x q[39];
x q[11];
x q[5];
x q[48];
x q[3];
x q[29];
x q[42];
x q[45];
x q[32];
x q[25];
x q[35];
x q[50];
x q[1];
cx q[18], q[13];
cx q[37], q[31];
cx q[6], q[12];
cx q[30], q[36];
cx q[20], q[26];
cx q[14], q[19];
cx q[34], q[28];
cx q[27], q[21];
cx q[2], q[7];
cx q[9], q[4];
cx q[40], q[46];
cx q[49], q[43];
cx q[47], q[41];
cx q[17], q[23];
x q[0];
x q[8];
x q[51];
x q[52];
x q[3];
x q[41];
x q[23];
x q[48];
x q[34];
x q[35];
x q[9];
x q[11];
x q[1];
x q[33];
x q[7];
x q[28];
x q[39];
x q[30];
x q[40];
x q[26];
x q[45];
x q[2];
x q[46];
x q[4];
x q[17];
x q[12];
x q[32];
x q[20];
x q[5];
cx q[44], q[38];
cx q[37], q[31];
cx q[53], q[47];
cx q[10], q[16];
cx q[43], q[50];
cx q[15], q[21];
cx q[14], q[19];
cx q[36], q[42];
cx q[25], q[18];
cx q[22], q[29];
x q[49];
x q[13];
cx q[6], q[0];
x q[8];
x q[51];
x q[52];
x q[38];
x q[46];
x q[14];
x q[15];
x q[30];
x q[1];
x q[32];
x q[19];
x q[11];
x q[25];
x q[4];
x q[9];
x q[17];
x q[12];
x q[50];
x q[3];
x q[47];
x q[37];
x q[23];
x q[5];
x q[2];
x q[53];
x q[7];
x q[42];
x q[43];
x q[36];
x q[20];
x q[41];
x q[48];
x q[44];
cx q[22], q[16];
cx q[33], q[28];
cx q[24], q[18];
cx q[27], q[21];
cx q[39], q[45];
cx q[40], q[34];
cx q[35], q[29];
cx q[31], q[26];
x q[0];
x q[13];
x q[10];
x q[6];
x q[52];
x q[26];
x q[11];
x q[2];
x q[12];
x q[36];
x q[31];
x q[27];
x q[53];
x q[23];
x q[50];
x q[21];
x q[42];
x q[4];
x q[45];
x q[48];
x q[9];
x q[19];
x q[5];
x q[14];
x q[24];
x q[37];
x q[17];
cx q[44], q[38];
cx q[25], q[30];
cx q[40], q[46];
cx q[33], q[39];
cx q[49], q[43];
cx q[15], q[20];
cx q[35], q[29];
cx q[47], q[41];
cx q[7], q[1];
cx q[22], q[28];
x q[10];
x q[3];
x q[0];
x q[13];
x q[18];
x q[52];
x q[29];
x q[42];
x q[36];
x q[48];
x q[31];
x q[38];
x q[12];
x q[22];
x q[21];
x q[47];
x q[5];
x q[23];
x q[9];
x q[37];
x q[1];
x q[15];
x q[7];
x q[4];
cx q[50], q[44];
cx q[2], q[8];
cx q[14], q[19];
cx q[34], q[28];
cx q[30], q[24];
cx q[32], q[26];
cx q[53], q[46];
cx q[40], q[45];
cx q[17], q[11];
cx q[35], q[41];
cx q[33], q[27];
cx q[49], q[43];
x q[3];
x q[25];
x q[18];
x q[0];
x q[39];
x q[13];
cx q[10], q[16];
x q[52];
x q[44];
x q[49];
x q[41];
x q[30];
x q[2];
x q[36];
x q[19];
x q[22];
x q[8];
x q[24];
x q[1];
x q[48];
x q[15];
x q[14];
x q[7];
x q[50];
x q[43];
x q[4];
x q[9];
x q[40];
x q[47];
x q[38];
x q[21];
cx q[53], q[46];
cx q[27], q[32];
cx q[11], q[5];
cx q[20], q[26];
cx q[34], q[28];
cx q[37], q[31];
cx q[6], q[12];
cx q[45], q[51];
cx q[35], q[29];
cx q[17], q[23];
x q[3];
x q[13];
x q[39];
x q[33];
x q[25];
x q[16];
x q[0];
x q[52];
x q[44];
x q[51];
x q[36];
x q[21];
x q[30];
x q[7];
x q[37];
x q[38];
x q[2];
x q[19];
x q[24];
x q[6];
x q[32];
x q[27];
x q[8];
x q[29];
x q[20];
x q[14];
x q[35];
cx q[46], q[41];
cx q[12], q[18];
cx q[11], q[5];
cx q[22], q[28];
cx q[15], q[9];
cx q[17], q[23];
cx q[49], q[43];
cx q[53], q[47];
cx q[31], q[26];
cx q[42], q[48];
cx q[40], q[34];
x q[25];
x q[1];
x q[33];
cx q[10], q[16];
x q[0];
x q[3];
x q[52];
x q[29];
x q[42];
x q[34];
x q[38];
x q[48];
x q[46];
x q[12];
x q[27];
x q[23];
x q[35];
x q[19];
x q[14];
x q[41];
x q[49];
x q[7];
x q[26];
x q[32];
x q[36];
cx q[39], q[44];
cx q[15], q[20];
cx q[22], q[17];
cx q[37], q[31];
cx q[28], q[21];
cx q[45], q[51];
cx q[18], q[13];
cx q[30], q[24];
cx q[43], q[50];
cx q[2], q[8];
cx q[53], q[47];
cx q[9], q[4];
x q[6];
x q[1];
x q[5];
x q[10];
x q[0];
x q[52];
x q[44];
x q[4];
x q[8];
x q[53];
x q[42];
x q[48];
x q[29];
x q[45];
x q[41];
x q[9];
x q[32];
x q[36];
x q[22];
x q[30];
x q[14];
x q[49];
x q[37];
x q[47];
x q[15];
x q[35];
x q[28];
x q[23];
x q[34];
x q[51];
x q[12];
cx q[17], q[11];
cx q[20], q[27];
cx q[43], q[50];
cx q[25], q[18];
cx q[2], q[7];
cx q[16], q[21];
cx q[40], q[46];
cx q[13], q[19];
cx q[31], q[38];
x q[5];
cx q[6], q[1];
cx q[33], q[39];
x q[0];
x q[26];
x q[44];
x q[4];
x q[23];
x q[12];
x q[35];
x q[46];
x q[7];
x q[25];
x q[31];
x q[53];
x q[51];
x q[17];
x q[45];
x q[2];
x q[18];
x q[36];
x q[15];
x q[13];
x q[16];
x q[34];
x q[9];
x q[20];
x q[29];
x q[48];
x q[27];
cx q[49], q[43];
cx q[37], q[42];
cx q[32], q[38];
cx q[14], q[8];
cx q[30], q[24];
cx q[22], q[28];
x q[19];
x q[40];
x q[39];
x q[6];
cx q[10], q[5];
cx q[47], q[41];
x q[1];
x q[0];
x q[44];
x q[7];
x q[35];
x q[51];
x q[18];
x q[34];
x q[27];
x q[30];
x q[12];
x q[13];
x q[9];
x q[20];
x q[45];
x q[24];
x q[4];
x q[14];
x q[36];
x q[15];
x q[8];
x q[37];
x q[32];
x q[49];
x q[2];
x q[48];
x q[38];
cx q[43], q[50];
cx q[17], q[23];
cx q[25], q[31];
cx q[22], q[29];
cx q[33], q[28];
cx q[53], q[46];
cx q[16], q[21];
x q[40];
x q[19];
x q[39];
x q[47];
cx q[10], q[5];
x q[1];
x q[0];
x q[44];
x q[31];
x q[50];
x q[43];
x q[48];
x q[7];
x q[51];
x q[53];
x q[4];
x q[49];
x q[37];
x q[23];
x q[8];
x q[29];
x q[2];
x q[13];
x q[35];
x q[45];
x q[28];
cx q[6], q[12];
cx q[46], q[41];
cx q[33], q[27];
cx q[17], q[11];
cx q[20], q[14];
cx q[22], q[16];
cx q[15], q[21];
cx q[24], q[18];
cx q[32], q[38];
cx q[25], q[30];
cx q[9], q[3];
cx q[36], q[42];
x q[39];
x q[34];
x q[47];
x q[0];
x q[5];
x q[1];
x q[44];
x q[12];
x q[11];
x q[49];
x q[7];
x q[6];
x q[3];
x q[9];
x q[23];
x q[41];
x q[46];
x q[17];
x q[35];
x q[32];
x q[18];
x q[30];
x q[24];
x q[48];
x q[29];
x q[2];
x q[36];
x q[14];
x q[51];
x q[38];
x q[53];
x q[25];
x q[20];
cx q[31], q[26];
cx q[13], q[19];
cx q[40], q[45];
cx q[43], q[50];
cx q[15], q[8];
cx q[10], q[4];
cx q[22], q[28];
cx q[33], q[27];
x q[37];
x q[34];
x q[39];
x q[0];
x q[5];
x q[1];
x q[44];
x q[9];
x q[19];
x q[43];
x q[36];
x q[38];
x q[17];
x q[18];
x q[30];
x q[4];
x q[8];
x q[6];
x q[28];
x q[10];
x q[12];
x q[3];
x q[24];
x q[48];
x q[35];
x q[7];
x q[2];
x q[53];
x q[50];
cx q[33], q[40];
cx q[45], q[52];
cx q[31], q[26];
cx q[47], q[41];
cx q[15], q[20];
cx q[49], q[42];
cx q[22], q[29];
x q[32];
x q[25];
x q[13];
x q[14];
x q[23];
x q[51];
x q[0];
x q[5];
x q[1];
x q[31];
x q[3];
x q[52];
x q[33];
x q[36];
x q[12];
x q[8];
x q[49];
x q[38];
x q[2];
x q[43];
x q[26];
x q[19];
x q[7];
cx q[50], q[44];
cx q[40], q[34];
cx q[20], q[27];
cx q[35], q[29];
cx q[15], q[9];
cx q[39], q[45];
cx q[24], q[18];
cx q[53], q[47];
cx q[22], q[16];
cx q[17], q[11];
cx q[37], q[30];
cx q[28], q[21];
cx q[46], q[41];
cx q[10], q[4];
cx q[42], q[48];
x q[14];
x q[25];
x q[23];
x q[6];
x q[13];
x q[32];
x q[0];
x q[5];
x q[1];
x q[50];
x q[22];
x q[34];
x q[49];
x q[48];
x q[43];
x q[18];
x q[4];
x q[40];
x q[12];
x q[41];
x q[3];
x q[35];
x q[10];
x q[24];
x q[47];
x q[53];
x q[8];
x q[52];
x q[21];
x q[28];
x q[33];
x q[39];
cx q[17], q[11];
cx q[15], q[20];
cx q[2], q[7];
cx q[37], q[30];
cx q[44], q[51];
cx q[9], q[16];
cx q[26], q[19];
cx q[36], q[42];
x q[38];
x q[32];
x q[6];
x q[45];
x q[31];
x q[0];
x q[1];
x q[10];
x q[21];
x q[36];
x q[4];
x q[49];
x q[44];
x q[51];
x q[39];
x q[47];
x q[8];
x q[2];
x q[33];
x q[16];
x q[12];
x q[37];
x q[3];
cx q[43], q[50];
cx q[42], q[48];
cx q[25], q[18];
cx q[30], q[24];
cx q[53], q[46];
cx q[35], q[41];
cx q[14], q[19];
cx q[40], q[34];
cx q[15], q[9];
cx q[20], q[26];
cx q[22], q[29];
cx q[7], q[13];
cx q[17], q[23];
x q[11];
x q[52];
x q[28];
x q[32];
cx q[31], q[38];
x q[0];
x q[1];
x q[18];
x q[41];
x q[29];
x q[48];
x q[8];
x q[35];
x q[36];
x q[34];
x q[30];
x q[26];
x q[20];
x q[24];
x q[53];
x q[37];
x q[15];
x q[16];
x q[13];
x q[46];
x q[19];
x q[14];
x q[10];
x q[4];
x q[25];
cx q[43], q[50];
cx q[17], q[23];
cx q[33], q[39];
cx q[9], q[3];
cx q[6], q[12];
cx q[2], q[7];
cx q[27], q[21];
cx q[45], q[51];
cx q[49], q[42];
x q[38];
x q[31];
x q[28];
x q[44];
x q[11];
x q[52];
x q[22];
x q[32];
x q[47];
x q[0];
x q[1];
x q[43];
x q[23];
x q[46];
x q[37];
x q[4];
x q[51];
x q[17];
x q[53];
x q[3];
x q[42];
x q[50];
x q[9];
x q[36];
x q[19];
x q[24];
x q[26];
x q[18];
cx q[35], q[29];
cx q[15], q[8];
cx q[25], q[30];
cx q[6], q[12];
cx q[34], q[41];
cx q[27], q[21];
cx q[7], q[13];
cx q[33], q[40];
cx q[39], q[45];
cx q[20], q[14];
x q[47];
x q[31];
x q[38];
x q[16];
x q[52];
cx q[10], q[5];
x q[48];
x q[32];
x q[11];
x q[0];
x q[1];
x q[43];
x q[18];
x q[21];
x q[24];
x q[37];
x q[51];
x q[36];
x q[30];
x q[3];
x q[45];
x q[8];
x q[23];
x q[14];
x q[12];
x q[27];
cx q[7], q[13];
cx q[25], q[19];
cx q[34], q[29];
cx q[9], q[4];
cx q[39], q[44];
cx q[20], q[26];
cx q[53], q[46];
cx q[35], q[41];
cx q[22], q[17];
cx q[49], q[42];
cx q[33], q[28];
x q[47];
x q[5];
x q[38];
x q[52];
x q[6];
x q[16];
x q[48];
x q[11];
x q[0];
x q[1];
x q[43];
x q[27];
x q[20];
x q[29];
x q[41];
x q[53];
x q[44];
x q[49];
x q[26];
x q[33];
x q[35];
x q[30];
x q[12];
x q[24];
x q[51];
cx q[17], q[23];
cx q[39], q[45];
cx q[15], q[21];
cx q[25], q[18];
cx q[37], q[31];
cx q[22], q[28];
cx q[36], q[42];
cx q[7], q[13];
cx q[9], q[3];
cx q[14], q[19];
cx q[40], q[46];
cx q[10], q[4];
cx q[2], q[8];
x q[34];
x q[5];
x q[6];
cx q[32], q[38];
x q[48];
x q[0];
x q[10];
x q[12];
x q[28];
x q[46];
x q[14];
x q[36];
x q[51];
x q[22];
x q[49];
x q[37];
x q[40];
x q[9];
x q[42];
x q[20];
x q[3];
x q[35];
x q[4];
x q[39];
x q[7];
x q[53];
x q[2];
x q[44];
x q[29];
x q[19];
x q[45];
cx q[43], q[50];
cx q[33], q[27];
cx q[15], q[21];
cx q[18], q[13];
cx q[17], q[23];
cx q[47], q[41];
cx q[25], q[31];
cx q[30], q[24];
x q[8];
x q[26];
x q[32];
x q[38];
x q[6];
x q[48];
x q[0];
x q[12];
x q[4];
x q[46];
x q[51];
x q[22];
x q[28];
x q[19];
x q[2];
x q[27];
x q[14];
x q[23];
x q[31];
x q[18];
x q[41];
x q[3];
x q[24];
x q[36];
cx q[50], q[44];
cx q[17], q[11];
cx q[35], q[29];
cx q[33], q[39];
cx q[10], q[5];
cx q[49], q[43];
cx q[25], q[30];
cx q[7], q[13];
cx q[16], q[21];
cx q[45], q[52];
cx q[15], q[9];
cx q[40], q[34];
cx q[53], q[47];
cx q[37], q[42];
x q[38];
x q[8];
x q[6];
x q[20];
x q[48];
x q[0];
x q[44];
x q[17];
x q[51];
x q[37];
x q[30];
x q[22];
x q[34];
x q[45];
x q[50];
x q[47];
x q[43];
x q[9];
x q[16];
x q[3];
x q[15];
x q[5];
x q[2];
x q[39];
x q[7];
x q[10];
x q[52];
x q[25];
x q[4];
x q[36];
x q[11];
x q[49];
x q[31];
x q[53];
x q[41];
x q[46];
x q[42];
x q[33];
x q[14];
x q[23];
x q[24];
cx q[12], q[18];
cx q[28], q[21];
cx q[27], q[32];
cx q[35], q[29];
cx q[26], q[19];
x q[13];
x q[40];
x q[0];
x q[46];
x q[23];
x q[18];
x q[35];
x q[14];
x q[36];
x q[41];
x q[47];
x q[42];
x q[22];
x q[37];
x q[53];
x q[7];
x q[52];
cx q[50], q[44];
cx q[34], q[29];
cx q[15], q[20];
cx q[49], q[43];
cx q[31], q[26];
cx q[9], q[4];
cx q[27], q[21];
cx q[25], q[19];
cx q[10], q[16];
cx q[6], q[12];
cx q[45], q[51];
cx q[30], q[24];
cx q[33], q[28];
cx q[17], q[11];
cx q[8], q[3];
cx q[32], q[38];
x q[13];
x q[40];
x q[2];
x q[5];
x q[50];
x q[44];
x q[37];
x q[3];
x q[33];
x q[16];
x q[21];
x q[36];
x q[53];
x q[10];
x q[12];
x q[11];
x q[45];
x q[28];
x q[35];
x q[34];
x q[7];
x q[51];
x q[30];
x q[17];
x q[22];
cx q[47], q[41];
cx q[9], q[4];
cx q[32], q[38];
cx q[46], q[52];
cx q[49], q[43];
cx q[14], q[8];
cx q[20], q[27];
cx q[6], q[1];
cx q[24], q[18];
cx q[29], q[23];
cx q[42], q[48];
cx q[25], q[31];
x q[13];
x q[19];
x q[26];
x q[40];
x q[2];
x q[15];
x q[50];
x q[51];
x q[43];
x q[10];
x q[38];
x q[24];
x q[35];
x q[21];
x q[12];
x q[28];
x q[48];
x q[4];
x q[8];
x q[31];
x q[34];
x q[49];
x q[29];
x q[6];
x q[33];
x q[44];
x q[27];
x q[32];
x q[18];
x q[52];
x q[37];
x q[25];
cx q[17], q[23];
cx q[7], q[1];
cx q[20], q[14];
cx q[9], q[3];
cx q[39], q[45];
cx q[11], q[5];
cx q[53], q[47];
cx q[46], q[41];
cx q[30], q[36];
x q[16];
x q[15];
x q[13];
x q[19];
x q[50];
x q[9];
x q[23];
x q[49];
x q[51];
x q[25];
x q[12];
x q[31];
x q[21];
x q[45];
x q[44];
x q[52];
x q[41];
x q[10];
x q[36];
x q[28];
x q[11];
x q[35];
x q[30];
x q[47];
x q[3];
x q[29];
x q[5];
x q[4];
x q[38];
x q[53];
x q[8];
x q[34];
cx q[24], q[18];
cx q[39], q[32];
cx q[42], q[48];
cx q[22], q[17];
cx q[6], q[1];
cx q[20], q[26];
cx q[37], q[43];
cx q[7], q[14];
x q[27];
x q[16];
x q[13];
x q[33];
x q[15];
cx q[40], q[46];
x q[19];
x q[3];
x q[17];
x q[39];
x q[1];
x q[25];
x q[23];
x q[21];
x q[37];
x q[30];
x q[20];
x q[35];
x q[49];
x q[45];
x q[48];
x q[43];
x q[7];
x q[51];
x q[6];
x q[47];
x q[29];
x q[12];
x q[14];
x q[4];
x q[42];
x q[36];
cx q[50], q[44];
cx q[22], q[28];
cx q[2], q[8];
cx q[34], q[41];
cx q[31], q[26];
cx q[24], q[18];
cx q[11], q[5];
cx q[32], q[38];
x q[52];
x q[16];
x q[15];
x q[53];
x q[9];
x q[46];
cx q[33], q[27];
x q[19];
x q[26];
x q[3];
x q[47];
x q[4];
x q[20];
x q[21];
x q[37];
x q[23];
x q[18];
x q[25];
x q[36];
x q[2];
x q[32];
x q[43];
x q[29];
x q[49];
x q[24];
x q[28];
x q[42];
x q[12];
x q[8];
x q[51];
x q[1];
x q[41];
cx q[44], q[38];
cx q[40], q[34];
cx q[7], q[13];
cx q[39], q[45];
cx q[6], q[0];
cx q[17], q[11];
x q[53];
x q[35];
x q[15];
x q[50];
x q[16];
x q[14];
x q[9];
x q[46];
x q[30];
cx q[33], q[27];
x q[52];
x q[48];
x q[22];
x q[40];
x q[3];
x q[17];
x q[4];
x q[18];
x q[23];
x q[25];
x q[43];
x q[6];
x q[7];
x q[51];
x q[45];
x q[13];
x q[36];
x q[1];
x q[0];
cx q[39], q[44];
cx q[49], q[42];
cx q[34], q[29];
cx q[20], q[26];
cx q[37], q[31];
cx q[47], q[41];
cx q[11], q[5];
cx q[2], q[8];
cx q[28], q[21];
cx q[32], q[38];
x q[24];
x q[30];
x q[33];
x q[27];
x q[14];
x q[35];
x q[9];
cx q[10], q[16];
x q[44];
x q[21];
x q[25];
x q[51];
x q[23];
x q[31];
x q[17];
x q[28];
x q[0];
x q[38];
x q[3];
x q[29];
x q[32];
x q[4];
x q[13];
x q[36];
x q[6];
x q[49];
x q[1];
cx q[46], q[41];
cx q[11], q[5];
cx q[37], q[42];
cx q[43], q[50];
cx q[15], q[20];
cx q[12], q[18];
cx q[2], q[7];
cx q[39], q[45];
cx q[40], q[34];
cx q[53], q[47];

// measurement
measure q[25]->c[0];
measure q[49]->c[1];
measure q[15]->c[2];
measure q[33]->c[3];
measure q[2]->c[4];
measure q[22]->c[5];
measure q[10]->c[6];
measure q[7]->c[7];
measure q[6]->c[8];
measure q[17]->c[9];
measure q[37]->c[10];
measure q[39]->c[11];
measure q[30]->c[12];
measure q[43]->c[13];
measure q[20]->c[14];
measure q[27]->c[15];
measure q[36]->c[16];
measure q[24]->c[17];
measure q[32]->c[18];
measure q[35]->c[19];
measure q[40]->c[20];
measure q[11]->c[21];
measure q[50]->c[22];
measure q[12]->c[23];
measure q[18]->c[24];
measure q[45]->c[25];
measure q[53]->c[26];
measure q[34]->c[27];
measure q[47]->c[28];
measure q[31]->c[29];
measure q[26]->c[30];
measure q[28]->c[31];
measure q[1]->c[32];
measure q[0]->c[33];
measure q[14]->c[34];
measure q[13]->c[35];
measure q[44]->c[36];
measure q[8]->c[37];
measure q[9]->c[38];
measure q[38]->c[39];
measure q[42]->c[40];
measure q[46]->c[41];
measure q[51]->c[42];
measure q[41]->c[43];
measure q[16]->c[44];
measure q[19]->c[45];
measure q[48]->c[46];
measure q[5]->c[47];
measure q[29]->c[48];
measure q[4]->c[49];
measure q[52]->c[50];
measure q[23]->c[51];
measure q[21]->c[52];
measure q[3]->c[53];
