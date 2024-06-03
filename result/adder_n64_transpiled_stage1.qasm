OPENQASM 2.0;
include "qelib1.inc";
qreg q[127];
creg c[127];
x q[90];
x q[95];
x q[97];
x q[61];
x q[58];
x q[53];
x q[40];
x q[44];
x q[46];
x q[65];
x q[66];
x q[49];
x q[55];
x q[70];
x q[88];
x q[83];
x q[62];
x q[81];
x q[72];
x q[67];
x q[84];
x q[87];
x q[106];
x q[103];
x q[101];
x q[100];
x q[118];
x q[75];
x q[77];
cx q[76], q[75];
cx q[90], q[94];
cx q[95], q[96];
cx q[97], q[98];
cx q[61], q[60];
cx q[58], q[59];
cx q[53], q[41];
cx q[40], q[39];
cx q[44], q[45];
cx q[46], q[47];
cx q[66], q[73];
cx q[49], q[50];
cx q[70], q[74];
cx q[88], q[89];
cx q[83], q[82];
cx q[62], q[63];
cx q[81], q[80];
cx q[84], q[85];
cx q[87], q[93];
cx q[106], q[105];
cx q[103], q[102];
cx q[100], q[110];
cx q[118], q[119];
swap q[64], q[65];
swap q[67], q[68];
cx q[64], q[54];
cx q[68], q[69];
cx q[76], q[77];
swap q[72], q[81];
swap q[49], q[55];
swap q[92], q[102];
swap q[61], q[62];
swap q[58], q[59];
swap q[88], q[89];
swap q[86], q[87];
swap q[83], q[84];
cx q[49], q[48];
cx q[101], q[102];
rz(pi/2) q[76];
swap q[80], q[81];
swap q[60], q[61];
swap q[54], q[64];
sx q[76];
swap q[49], q[50];
swap q[92], q[102];
rz(pi/2) q[76];
cx q[76], q[75];
rz(-pi/4) q[76];
cx q[76], q[77];
rz(pi/4) q[76];
cx q[76], q[75];
rz(-pi/4) q[76];
rz(pi/4) q[75];
cx q[76], q[77];
rz(3*pi/4) q[76];
sx q[76];
rz(pi/2) q[76];
swap q[75], q[76];
cx q[75], q[90];
cx q[76], q[77];
rz(pi/2) q[90];
rz(-pi/4) q[76];
rz(pi/4) q[77];
sx q[90];
cx q[76], q[77];
rz(pi/2) q[90];
cx q[90], q[94];
rz(-pi/4) q[90];
cx q[75], q[90];
rz(pi/4) q[90];
cx q[90], q[94];
rz(-pi/4) q[90];
rz(pi/4) q[94];
cx q[75], q[90];
rz(3*pi/4) q[90];
sx q[90];
rz(pi/2) q[90];
swap q[90], q[94];
cx q[94], q[95];
cx q[75], q[90];
rz(pi/2) q[95];
rz(pi/4) q[75];
rz(-pi/4) q[90];
sx q[95];
cx q[75], q[90];
rz(pi/2) q[95];
cx q[95], q[96];
rz(-pi/4) q[95];
cx q[94], q[95];
rz(pi/4) q[95];
cx q[95], q[96];
rz(-pi/4) q[95];
rz(pi/4) q[96];
cx q[94], q[95];
rz(3*pi/4) q[95];
sx q[95];
rz(pi/2) q[95];
swap q[95], q[96];
cx q[96], q[97];
cx q[94], q[95];
rz(pi/2) q[97];
rz(pi/4) q[94];
rz(-pi/4) q[95];
sx q[97];
cx q[94], q[95];
rz(pi/2) q[97];
cx q[97], q[98];
rz(-pi/4) q[97];
cx q[96], q[97];
rz(pi/4) q[97];
cx q[97], q[98];
rz(-pi/4) q[97];
rz(pi/4) q[98];
cx q[96], q[97];
rz(3*pi/4) q[97];
sx q[97];
rz(pi/2) q[97];
swap q[97], q[98];
cx q[96], q[97];
swap q[91], q[98];
rz(pi/4) q[96];
rz(-pi/4) q[97];
cx q[91], q[79];
cx q[96], q[97];
rz(pi/2) q[91];
swap q[79], q[80];
cx q[79], q[78];
sx q[91];
swap q[80], q[81];
rz(pi/2) q[91];
swap q[72], q[81];
cx q[62], q[72];
swap q[91], q[98];
cx q[98], q[97];
rz(pi/2) q[62];
rz(-pi/4) q[98];
sx q[62];
rz(pi/2) q[62];
swap q[97], q[98];
cx q[96], q[97];
cx q[62], q[61];
rz(pi/4) q[97];
rz(-pi/4) q[62];
cx q[97], q[98];
cx q[62], q[72];
rz(-pi/4) q[97];
rz(pi/4) q[98];
rz(pi/4) q[62];
cx q[96], q[97];
cx q[62], q[61];
rz(3*pi/4) q[97];
rz(pi/4) q[61];
rz(-pi/4) q[62];
sx q[97];
cx q[62], q[72];
rz(pi/2) q[97];
rz(3*pi/4) q[62];
sx q[62];
swap q[97], q[98];
cx q[96], q[97];
rz(pi/2) q[62];
rz(pi/4) q[96];
rz(-pi/4) q[97];
swap q[61], q[62];
cx q[96], q[97];
cx q[62], q[72];
swap q[60], q[61];
cx q[60], q[59];
rz(-pi/4) q[62];
rz(pi/4) q[72];
swap q[96], q[97];
cx q[97], q[98];
rz(pi/2) q[59];
cx q[62], q[72];
cx q[97], q[96];
sx q[59];
swap q[61], q[62];
rz(pi/2) q[97];
rz(pi/2) q[59];
sx q[97];
cx q[59], q[58];
rz(pi/2) q[97];
rz(-pi/4) q[59];
cx q[60], q[59];
swap q[96], q[97];
cx q[96], q[95];
rz(pi/4) q[59];
rz(-pi/4) q[96];
cx q[59], q[58];
rz(pi/4) q[58];
rz(-pi/4) q[59];
swap q[95], q[96];
cx q[94], q[95];
cx q[60], q[59];
rz(pi/4) q[95];
rz(3*pi/4) q[59];
cx q[95], q[96];
sx q[59];
rz(-pi/4) q[95];
rz(pi/4) q[96];
rz(pi/2) q[59];
cx q[94], q[95];
swap q[59], q[60];
rz(3*pi/4) q[95];
cx q[59], q[58];
cx q[60], q[53];
sx q[95];
rz(-pi/4) q[58];
rz(pi/4) q[59];
rz(pi/2) q[53];
rz(pi/2) q[95];
cx q[59], q[58];
sx q[53];
rz(pi/2) q[53];
swap q[95], q[96];
swap q[58], q[59];
cx q[94], q[95];
cx q[53], q[41];
rz(pi/4) q[94];
rz(-pi/4) q[95];
rz(-pi/4) q[53];
cx q[94], q[95];
cx q[60], q[53];
rz(pi/4) q[53];
swap q[94], q[95];
cx q[95], q[96];
cx q[53], q[41];
cx q[95], q[94];
rz(pi/4) q[41];
rz(-pi/4) q[53];
rz(pi/2) q[95];
cx q[60], q[53];
sx q[95];
rz(3*pi/4) q[53];
rz(pi/2) q[95];
sx q[53];
rz(pi/2) q[53];
swap q[94], q[95];
cx q[94], q[90];
swap q[41], q[53];
rz(-pi/4) q[94];
cx q[60], q[53];
cx q[41], q[40];
rz(-pi/4) q[53];
rz(pi/4) q[60];
rz(pi/2) q[40];
swap q[90], q[94];
cx q[75], q[90];
cx q[60], q[53];
sx q[40];
rz(pi/4) q[90];
rz(pi/2) q[40];
cx q[90], q[94];
cx q[40], q[39];
rz(-pi/4) q[90];
rz(pi/4) q[94];
rz(-pi/4) q[40];
cx q[75], q[90];
cx q[41], q[40];
rz(3*pi/4) q[90];
rz(pi/4) q[40];
sx q[90];
cx q[40], q[39];
rz(pi/2) q[90];
rz(pi/4) q[39];
rz(-pi/4) q[40];
cx q[41], q[40];
swap q[90], q[94];
cx q[75], q[90];
rz(3*pi/4) q[40];
rz(pi/4) q[75];
rz(-pi/4) q[90];
sx q[40];
cx q[75], q[90];
rz(pi/2) q[40];
swap q[75], q[90];
swap q[40], q[41];
cx q[90], q[94];
cx q[40], q[39];
cx q[41], q[42];
cx q[90], q[75];
rz(-pi/4) q[39];
rz(pi/4) q[40];
rz(pi/2) q[41];
swap q[42], q[43];
rz(pi/2) q[90];
cx q[40], q[39];
sx q[41];
cx q[44], q[43];
sx q[90];
rz(pi/2) q[41];
rz(pi/2) q[44];
rz(pi/2) q[90];
sx q[44];
swap q[40], q[41];
cx q[40], q[39];
rz(pi/2) q[44];
swap q[75], q[90];
cx q[75], q[76];
rz(-pi/4) q[40];
cx q[44], q[45];
rz(-pi/4) q[75];
cx q[41], q[40];
rz(-pi/4) q[44];
rz(pi/4) q[40];
cx q[44], q[43];
swap q[75], q[76];
cx q[76], q[77];
cx q[40], q[39];
rz(pi/4) q[44];
rz(pi/4) q[76];
rz(pi/4) q[39];
rz(-pi/4) q[40];
cx q[44], q[45];
cx q[76], q[75];
cx q[41], q[40];
rz(pi/4) q[45];
rz(-pi/4) q[44];
rz(-pi/4) q[76];
rz(pi/4) q[75];
rz(3*pi/4) q[40];
cx q[44], q[43];
cx q[76], q[77];
sx q[40];
rz(3*pi/4) q[44];
rz(3*pi/4) q[76];
rz(pi/2) q[40];
sx q[44];
sx q[76];
rz(pi/2) q[44];
swap q[40], q[41];
rz(pi/2) q[76];
cx q[40], q[39];
swap q[44], q[45];
rz(-pi/4) q[39];
rz(pi/4) q[40];
cx q[44], q[43];
cx q[45], q[46];
swap q[76], q[77];
cx q[75], q[76];
cx q[40], q[39];
rz(-pi/4) q[44];
rz(pi/4) q[43];
rz(pi/2) q[46];
rz(-pi/4) q[75];
rz(pi/4) q[76];
cx q[40], q[41];
cx q[44], q[43];
sx q[46];
cx q[75], q[76];
cx q[40], q[39];
rz(pi/2) q[46];
cx q[77], q[76];
rz(pi/2) q[40];
cx q[46], q[47];
cx q[75], q[76];
sx q[40];
rz(-pi/4) q[46];
rz(pi/2) q[40];
cx q[45], q[46];
rz(pi/4) q[46];
swap q[40], q[41];
cx q[41], q[53];
cx q[46], q[47];
rz(-pi/4) q[41];
rz(pi/4) q[47];
rz(-pi/4) q[46];
cx q[45], q[46];
swap q[41], q[53];
cx q[60], q[53];
rz(3*pi/4) q[46];
rz(pi/4) q[53];
sx q[46];
cx q[53], q[41];
rz(pi/2) q[46];
rz(pi/4) q[41];
rz(-pi/4) q[53];
swap q[45], q[46];
cx q[60], q[53];
cx q[46], q[47];
cx q[45], q[54];
rz(3*pi/4) q[53];
rz(-pi/4) q[47];
rz(pi/4) q[46];
rz(pi/2) q[54];
sx q[53];
cx q[46], q[47];
sx q[54];
rz(pi/2) q[53];
rz(pi/2) q[54];
cx q[54], q[64];
swap q[53], q[60];
cx q[53], q[41];
rz(-pi/4) q[54];
rz(-pi/4) q[41];
rz(pi/4) q[53];
cx q[45], q[54];
cx q[53], q[41];
rz(pi/4) q[54];
cx q[53], q[60];
cx q[54], q[64];
cx q[53], q[41];
rz(-pi/4) q[54];
rz(pi/4) q[64];
rz(pi/2) q[53];
cx q[45], q[54];
sx q[53];
rz(3*pi/4) q[54];
rz(pi/2) q[53];
sx q[54];
rz(pi/2) q[54];
swap q[53], q[60];
cx q[60], q[59];
swap q[54], q[64];
rz(-pi/4) q[60];
cx q[45], q[54];
swap q[64], q[65];
cx q[65], q[66];
rz(-pi/4) q[54];
rz(pi/4) q[45];
swap q[59], q[60];
cx q[58], q[59];
rz(pi/2) q[66];
cx q[45], q[54];
rz(pi/4) q[59];
sx q[66];
cx q[59], q[60];
rz(pi/2) q[66];
rz(pi/4) q[60];
rz(-pi/4) q[59];
cx q[66], q[73];
cx q[58], q[59];
rz(-pi/4) q[66];
rz(3*pi/4) q[59];
cx q[65], q[66];
sx q[59];
rz(pi/4) q[66];
rz(pi/2) q[59];
cx q[66], q[73];
rz(-pi/4) q[66];
rz(pi/4) q[73];
swap q[58], q[59];
cx q[59], q[60];
cx q[65], q[66];
rz(-pi/4) q[60];
rz(pi/4) q[59];
rz(3*pi/4) q[66];
cx q[59], q[60];
sx q[66];
cx q[59], q[58];
rz(pi/2) q[66];
cx q[59], q[60];
cx q[66], q[67];
rz(pi/2) q[59];
rz(pi/2) q[66];
swap q[67], q[68];
sx q[59];
sx q[66];
cx q[55], q[68];
rz(pi/2) q[59];
rz(pi/2) q[66];
rz(pi/2) q[55];
sx q[55];
swap q[66], q[73];
swap q[59], q[60];
cx q[60], q[61];
cx q[65], q[66];
rz(pi/2) q[55];
rz(-pi/4) q[60];
rz(pi/4) q[65];
rz(-pi/4) q[66];
cx q[55], q[49];
cx q[65], q[66];
rz(-pi/4) q[55];
swap q[60], q[61];
cx q[73], q[66];
cx q[55], q[68];
rz(-pi/4) q[73];
rz(pi/4) q[55];
cx q[55], q[49];
swap q[66], q[73];
cx q[65], q[66];
rz(-pi/4) q[55];
rz(pi/4) q[49];
rz(pi/4) q[66];
cx q[55], q[68];
swap q[49], q[50];
cx q[66], q[73];
rz(3*pi/4) q[55];
rz(-pi/4) q[66];
rz(pi/4) q[73];
sx q[55];
cx q[65], q[66];
rz(pi/2) q[55];
rz(3*pi/4) q[66];
cx q[55], q[49];
sx q[66];
rz(pi/2) q[49];
rz(pi/2) q[66];
sx q[49];
rz(pi/2) q[49];
swap q[65], q[66];
cx q[66], q[73];
cx q[49], q[48];
rz(pi/4) q[66];
rz(-pi/4) q[73];
rz(-pi/4) q[49];
cx q[66], q[73];
cx q[55], q[49];
cx q[66], q[65];
rz(pi/4) q[49];
cx q[66], q[73];
cx q[49], q[48];
rz(pi/2) q[66];
rz(-pi/4) q[49];
rz(pi/4) q[48];
sx q[66];
cx q[55], q[49];
rz(pi/2) q[66];
rz(3*pi/4) q[49];
sx q[49];
swap q[65], q[66];
rz(pi/2) q[49];
swap q[64], q[65];
swap q[66], q[67];
cx q[64], q[54];
swap q[49], q[55];
rz(-pi/4) q[64];
cx q[49], q[48];
swap q[55], q[68];
swap q[45], q[54];
rz(pi/4) q[49];
rz(-pi/4) q[48];
cx q[54], q[64];
swap q[68], q[69];
cx q[69], q[70];
cx q[49], q[48];
rz(pi/4) q[64];
swap q[67], q[68];
rz(pi/2) q[70];
swap q[49], q[55];
swap q[54], q[64];
sx q[70];
cx q[50], q[49];
cx q[54], q[45];
rz(pi/2) q[70];
rz(-pi/4) q[50];
rz(pi/4) q[49];
rz(-pi/4) q[54];
rz(pi/4) q[45];
cx q[70], q[74];
cx q[50], q[49];
cx q[64], q[54];
rz(-pi/4) q[70];
rz(3*pi/4) q[54];
swap q[48], q[49];
cx q[69], q[70];
sx q[54];
swap q[49], q[55];
rz(pi/4) q[70];
rz(pi/2) q[54];
cx q[70], q[74];
swap q[54], q[64];
rz(-pi/4) q[70];
rz(pi/4) q[74];
cx q[54], q[45];
cx q[69], q[70];
rz(-pi/4) q[45];
rz(pi/4) q[54];
rz(3*pi/4) q[70];
cx q[54], q[45];
sx q[70];
cx q[54], q[64];
rz(pi/2) q[70];
cx q[54], q[45];
rz(pi/2) q[54];
swap q[70], q[74];
cx q[74], q[89];
cx q[69], q[70];
sx q[54];
rz(pi/2) q[89];
rz(pi/4) q[69];
rz(-pi/4) q[70];
rz(pi/2) q[54];
sx q[89];
cx q[69], q[70];
swap q[45], q[54];
rz(pi/2) q[89];
swap q[45], q[46];
cx q[89], q[88];
cx q[46], q[47];
rz(-pi/4) q[89];
rz(-pi/4) q[46];
cx q[74], q[89];
cx q[45], q[46];
rz(pi/4) q[89];
rz(pi/4) q[46];
cx q[89], q[88];
cx q[46], q[47];
rz(-pi/4) q[89];
rz(pi/4) q[88];
rz(pi/4) q[47];
rz(-pi/4) q[46];
cx q[74], q[89];
cx q[45], q[46];
rz(3*pi/4) q[89];
rz(3*pi/4) q[46];
sx q[89];
sx q[46];
rz(pi/2) q[89];
rz(pi/2) q[46];
swap q[88], q[89];
swap q[46], q[47];
cx q[74], q[89];
cx q[45], q[46];
cx q[88], q[87];
rz(pi/4) q[74];
rz(-pi/4) q[89];
rz(-pi/4) q[46];
rz(pi/4) q[45];
rz(pi/2) q[88];
swap q[86], q[87];
cx q[74], q[89];
cx q[45], q[46];
sx q[88];
swap q[85], q[86];
rz(pi/2) q[88];
cx q[84], q[85];
swap q[45], q[46];
cx q[46], q[47];
rz(pi/2) q[84];
cx q[88], q[89];
cx q[46], q[45];
sx q[84];
rz(-pi/4) q[88];
rz(pi/2) q[46];
rz(pi/2) q[84];
swap q[88], q[89];
sx q[46];
cx q[74], q[89];
swap q[83], q[84];
rz(pi/2) q[46];
rz(pi/4) q[89];
cx q[83], q[82];
swap q[84], q[85];
cx q[89], q[88];
rz(-pi/4) q[83];
swap q[45], q[46];
cx q[45], q[44];
rz(-pi/4) q[89];
rz(pi/4) q[88];
cx q[83], q[84];
rz(-pi/4) q[45];
cx q[74], q[89];
rz(pi/4) q[83];
rz(3*pi/4) q[89];
cx q[83], q[82];
swap q[44], q[45];
cx q[44], q[43];
sx q[89];
rz(-pi/4) q[83];
rz(pi/4) q[82];
rz(pi/4) q[44];
rz(pi/2) q[89];
cx q[83], q[84];
cx q[44], q[45];
rz(3*pi/4) q[83];
swap q[74], q[89];
rz(pi/4) q[45];
rz(-pi/4) q[44];
cx q[89], q[88];
sx q[83];
cx q[44], q[43];
rz(pi/4) q[89];
rz(-pi/4) q[88];
rz(pi/2) q[83];
rz(3*pi/4) q[44];
cx q[89], q[88];
swap q[82], q[83];
sx q[44];
cx q[89], q[74];
cx q[83], q[84];
swap q[81], q[82];
rz(pi/2) q[44];
cx q[89], q[88];
rz(-pi/4) q[83];
rz(pi/4) q[84];
swap q[72], q[81];
rz(pi/2) q[89];
cx q[72], q[62];
cx q[83], q[84];
swap q[43], q[44];
swap q[81], q[82];
cx q[45], q[44];
sx q[89];
rz(pi/2) q[62];
rz(-pi/4) q[45];
rz(pi/4) q[44];
rz(pi/2) q[89];
sx q[62];
cx q[45], q[44];
rz(pi/2) q[62];
swap q[74], q[89];
cx q[43], q[44];
cx q[74], q[70];
cx q[62], q[63];
cx q[45], q[44];
rz(-pi/4) q[74];
rz(-pi/4) q[62];
swap q[69], q[70];
cx q[70], q[74];
cx q[72], q[62];
rz(pi/4) q[74];
rz(pi/4) q[62];
cx q[62], q[63];
swap q[70], q[74];
cx q[70], q[69];
rz(-pi/4) q[62];
rz(pi/4) q[63];
rz(-pi/4) q[70];
rz(pi/4) q[69];
cx q[72], q[62];
cx q[74], q[70];
rz(3*pi/4) q[62];
rz(3*pi/4) q[70];
sx q[62];
sx q[70];
rz(pi/2) q[62];
rz(pi/2) q[70];
swap q[62], q[72];
cx q[72], q[81];
cx q[62], q[63];
swap q[70], q[74];
cx q[70], q[69];
rz(pi/2) q[81];
rz(pi/4) q[62];
rz(-pi/4) q[63];
rz(pi/4) q[70];
rz(-pi/4) q[69];
sx q[81];
cx q[62], q[63];
cx q[70], q[69];
rz(pi/2) q[81];
cx q[70], q[74];
cx q[81], q[80];
cx q[70], q[69];
rz(-pi/4) q[81];
rz(pi/2) q[70];
cx q[72], q[81];
sx q[70];
rz(pi/4) q[81];
rz(pi/2) q[70];
cx q[81], q[80];
rz(-pi/4) q[81];
rz(pi/4) q[80];
swap q[69], q[70];
cx q[72], q[81];
swap q[68], q[69];
cx q[68], q[55];
rz(3*pi/4) q[81];
rz(-pi/4) q[68];
sx q[81];
rz(pi/2) q[81];
swap q[55], q[68];
cx q[49], q[55];
swap q[80], q[81];
rz(pi/4) q[55];
cx q[80], q[79];
cx q[72], q[81];
cx q[55], q[68];
rz(pi/2) q[79];
rz(pi/4) q[72];
rz(-pi/4) q[81];
rz(-pi/4) q[55];
rz(pi/4) q[68];
sx q[79];
cx q[72], q[81];
cx q[49], q[55];
rz(pi/2) q[79];
rz(3*pi/4) q[55];
cx q[79], q[78];
sx q[55];
rz(-pi/4) q[79];
rz(pi/2) q[55];
cx q[80], q[79];
rz(pi/4) q[79];
swap q[55], q[68];
cx q[49], q[55];
cx q[79], q[78];
rz(pi/4) q[49];
rz(-pi/4) q[55];
rz(-pi/4) q[79];
rz(pi/4) q[78];
cx q[49], q[55];
cx q[80], q[79];
rz(3*pi/4) q[79];
swap q[49], q[55];
cx q[55], q[68];
sx q[79];
cx q[55], q[49];
rz(pi/2) q[79];
rz(pi/2) q[55];
swap q[79], q[80];
sx q[55];
cx q[79], q[78];
swap q[80], q[81];
rz(pi/2) q[55];
rz(pi/4) q[79];
rz(-pi/4) q[78];
swap q[72], q[81];
cx q[79], q[78];
swap q[62], q[72];
swap q[81], q[82];
swap q[49], q[55];
cx q[49], q[50];
swap q[62], q[63];
swap q[72], q[81];
swap q[78], q[79];
rz(-pi/4) q[49];
swap q[63], q[64];
swap q[62], q[72];
cx q[61], q[62];
cx q[49], q[48];
cx q[64], q[65];
rz(pi/4) q[61];
rz(pi/4) q[49];
rz(pi/2) q[64];
cx q[66], q[65];
cx q[61], q[60];
cx q[49], q[50];
sx q[64];
rz(pi/2) q[66];
rz(pi/4) q[60];
rz(-pi/4) q[61];
rz(-pi/4) q[49];
rz(pi/4) q[50];
rz(pi/2) q[64];
sx q[66];
cx q[61], q[62];
cx q[49], q[48];
rz(pi/2) q[66];
swap q[63], q[64];
rz(3*pi/4) q[61];
rz(3*pi/4) q[49];
cx q[66], q[67];
sx q[61];
sx q[49];
rz(-pi/4) q[66];
rz(pi/2) q[61];
rz(pi/2) q[49];
cx q[66], q[65];
rz(pi/4) q[66];
swap q[61], q[62];
swap q[48], q[49];
cx q[60], q[61];
cx q[50], q[49];
cx q[66], q[67];
rz(-pi/4) q[60];
rz(pi/4) q[61];
rz(-pi/4) q[50];
rz(pi/4) q[49];
rz(-pi/4) q[66];
rz(pi/4) q[67];
cx q[60], q[61];
cx q[50], q[49];
cx q[66], q[65];
cx q[62], q[61];
cx q[48], q[49];
rz(3*pi/4) q[66];
cx q[60], q[61];
cx q[50], q[49];
sx q[66];
swap q[62], q[63];
rz(pi/2) q[66];
swap q[62], q[72];
swap q[66], q[73];
swap q[72], q[81];
cx q[73], q[85];
swap q[66], q[67];
swap q[80], q[81];
cx q[80], q[79];
rz(pi/2) q[85];
cx q[66], q[65];
rz(-pi/4) q[80];
sx q[85];
rz(-pi/4) q[66];
rz(pi/4) q[65];
rz(pi/2) q[85];
cx q[66], q[65];
swap q[79], q[80];
cx q[78], q[79];
cx q[85], q[86];
rz(pi/4) q[79];
rz(-pi/4) q[85];
cx q[79], q[80];
cx q[73], q[85];
rz(-pi/4) q[79];
rz(pi/4) q[80];
rz(pi/4) q[85];
cx q[78], q[79];
cx q[85], q[86];
rz(3*pi/4) q[79];
rz(-pi/4) q[85];
rz(pi/4) q[86];
sx q[79];
cx q[73], q[85];
swap q[86], q[87];
rz(pi/2) q[79];
rz(3*pi/4) q[85];
sx q[85];
swap q[78], q[79];
cx q[79], q[80];
rz(pi/2) q[85];
rz(pi/4) q[79];
rz(-pi/4) q[80];
cx q[85], q[86];
cx q[79], q[80];
rz(pi/2) q[86];
cx q[79], q[78];
sx q[86];
cx q[79], q[80];
rz(pi/2) q[86];
rz(pi/2) q[79];
swap q[86], q[87];
sx q[79];
cx q[87], q[93];
swap q[85], q[86];
rz(pi/2) q[79];
cx q[73], q[85];
rz(-pi/4) q[87];
rz(pi/4) q[73];
rz(-pi/4) q[85];
cx q[86], q[87];
swap q[79], q[80];
cx q[80], q[81];
cx q[73], q[85];
rz(pi/4) q[87];
rz(-pi/4) q[80];
cx q[87], q[93];
rz(-pi/4) q[87];
rz(pi/4) q[93];
swap q[80], q[81];
cx q[82], q[81];
cx q[86], q[87];
rz(pi/4) q[81];
rz(3*pi/4) q[87];
cx q[81], q[80];
sx q[87];
rz(-pi/4) q[81];
rz(pi/4) q[80];
rz(pi/2) q[87];
cx q[82], q[81];
swap q[87], q[93];
rz(3*pi/4) q[81];
cx q[93], q[106];
cx q[86], q[87];
sx q[81];
rz(pi/2) q[106];
rz(pi/4) q[86];
rz(-pi/4) q[87];
rz(pi/2) q[81];
sx q[106];
cx q[86], q[87];
rz(pi/2) q[106];
swap q[81], q[82];
cx q[81], q[80];
cx q[106], q[105];
rz(pi/4) q[81];
rz(-pi/4) q[80];
rz(-pi/4) q[106];
cx q[81], q[80];
cx q[93], q[106];
cx q[81], q[82];
rz(pi/4) q[106];
cx q[81], q[80];
cx q[106], q[105];
rz(pi/2) q[81];
rz(-pi/4) q[106];
rz(pi/4) q[105];
sx q[81];
cx q[93], q[106];
rz(pi/2) q[81];
rz(3*pi/4) q[106];
sx q[106];
swap q[72], q[81];
cx q[72], q[62];
rz(pi/2) q[106];
rz(-pi/4) q[72];
swap q[105], q[106];
cx q[81], q[72];
cx q[93], q[106];
cx q[105], q[104];
rz(pi/4) q[72];
rz(pi/4) q[93];
rz(-pi/4) q[106];
rz(pi/2) q[105];
cx q[103], q[104];
cx q[72], q[62];
cx q[93], q[106];
sx q[105];
rz(pi/2) q[103];
rz(-pi/4) q[72];
rz(pi/4) q[62];
rz(pi/2) q[105];
sx q[103];
cx q[81], q[72];
rz(pi/2) q[103];
cx q[105], q[106];
rz(3*pi/4) q[72];
rz(-pi/4) q[105];
cx q[103], q[102];
sx q[72];
rz(-pi/4) q[103];
swap q[105], q[106];
rz(pi/2) q[72];
cx q[93], q[106];
cx q[103], q[104];
rz(pi/4) q[106];
rz(pi/4) q[103];
swap q[72], q[81];
cx q[72], q[62];
cx q[106], q[105];
cx q[103], q[102];
rz(pi/4) q[72];
rz(-pi/4) q[62];
rz(-pi/4) q[106];
rz(pi/4) q[105];
rz(-pi/4) q[103];
rz(pi/4) q[102];
cx q[72], q[62];
cx q[93], q[106];
cx q[103], q[104];
swap q[101], q[102];
cx q[72], q[81];
rz(3*pi/4) q[106];
rz(3*pi/4) q[103];
cx q[72], q[62];
sx q[106];
sx q[103];
rz(pi/2) q[72];
rz(pi/2) q[106];
rz(pi/2) q[103];
sx q[72];
cx q[103], q[102];
swap q[93], q[106];
rz(pi/2) q[72];
cx q[106], q[105];
rz(pi/2) q[102];
rz(pi/4) q[106];
rz(-pi/4) q[105];
sx q[102];
swap q[72], q[81];
cx q[106], q[105];
rz(pi/2) q[102];
swap q[81], q[82];
cx q[82], q[83];
cx q[106], q[93];
cx q[102], q[92];
rz(-pi/4) q[82];
cx q[106], q[105];
rz(-pi/4) q[102];
rz(pi/2) q[106];
cx q[103], q[102];
swap q[82], q[83];
cx q[83], q[84];
sx q[106];
rz(pi/4) q[102];
rz(pi/4) q[83];
rz(pi/2) q[106];
cx q[102], q[92];
cx q[83], q[82];
rz(-pi/4) q[102];
rz(pi/4) q[92];
swap q[93], q[106];
rz(-pi/4) q[83];
rz(pi/4) q[82];
cx q[93], q[87];
cx q[103], q[102];
cx q[83], q[84];
rz(-pi/4) q[93];
rz(3*pi/4) q[102];
rz(3*pi/4) q[83];
sx q[102];
swap q[87], q[93];
sx q[83];
cx q[86], q[87];
rz(pi/2) q[102];
rz(pi/2) q[83];
rz(pi/4) q[87];
swap q[101], q[102];
cx q[87], q[93];
cx q[101], q[100];
swap q[102], q[103];
swap q[83], q[84];
cx q[82], q[83];
rz(-pi/4) q[87];
rz(pi/4) q[93];
rz(pi/2) q[100];
cx q[102], q[92];
cx q[103], q[104];
rz(-pi/4) q[82];
rz(pi/4) q[83];
cx q[86], q[87];
sx q[100];
rz(pi/4) q[102];
rz(-pi/4) q[92];
rz(-pi/4) q[103];
rz(pi/4) q[104];
cx q[82], q[83];
rz(3*pi/4) q[87];
rz(pi/2) q[100];
cx q[102], q[92];
cx q[103], q[104];
cx q[84], q[83];
sx q[87];
cx q[100], q[110];
swap q[92], q[102];
cx q[82], q[83];
rz(pi/2) q[87];
rz(-pi/4) q[100];
cx q[101], q[100];
swap q[86], q[87];
cx q[87], q[93];
rz(pi/4) q[100];
rz(pi/4) q[87];
rz(-pi/4) q[93];
cx q[100], q[110];
cx q[87], q[93];
rz(-pi/4) q[100];
rz(pi/4) q[110];
cx q[87], q[86];
cx q[101], q[100];
cx q[87], q[93];
rz(3*pi/4) q[100];
rz(pi/2) q[87];
sx q[100];
sx q[87];
rz(pi/2) q[100];
rz(pi/2) q[87];
swap q[100], q[110];
cx q[110], q[118];
cx q[101], q[100];
swap q[86], q[87];
cx q[86], q[85];
rz(pi/2) q[118];
rz(pi/4) q[101];
rz(-pi/4) q[100];
rz(-pi/4) q[86];
sx q[118];
cx q[101], q[100];
rz(pi/2) q[118];
swap q[85], q[86];
cx q[73], q[85];
cx q[118], q[119];
rz(pi/4) q[85];
rz(-pi/4) q[118];
cx q[85], q[86];
cx q[110], q[118];
rz(-pi/4) q[85];
rz(pi/4) q[86];
rz(pi/4) q[118];
cx q[73], q[85];
cx q[118], q[119];
rz(3*pi/4) q[85];
rz(-pi/4) q[118];
rz(pi/4) q[119];
sx q[85];
cx q[110], q[118];
rz(pi/2) q[85];
rz(3*pi/4) q[118];
sx q[118];
swap q[73], q[85];
cx q[85], q[86];
rz(pi/2) q[118];
rz(pi/4) q[85];
rz(-pi/4) q[86];
cx q[118], q[117];
cx q[85], q[86];
rz(pi/2) q[118];
cx q[85], q[73];
sx q[118];
cx q[85], q[86];
rz(pi/2) q[118];
rz(pi/2) q[85];
swap q[110], q[118];
sx q[85];
cx q[118], q[119];
rz(pi/2) q[85];
rz(pi/4) q[118];
rz(-pi/4) q[119];
cx q[118], q[119];
swap q[73], q[85];
cx q[73], q[66];
swap q[110], q[118];
rz(-pi/4) q[73];
cx q[118], q[119];
rz(-pi/4) q[118];
swap q[66], q[73];
cx q[66], q[65];
cx q[110], q[118];
rz(pi/4) q[66];
rz(pi/4) q[118];
cx q[66], q[73];
cx q[118], q[119];
rz(-pi/4) q[66];
rz(pi/4) q[73];
rz(-pi/4) q[118];
rz(pi/4) q[119];
cx q[66], q[65];
cx q[110], q[118];
rz(3*pi/4) q[66];
rz(3*pi/4) q[118];
sx q[66];
sx q[118];
rz(pi/2) q[66];
rz(pi/2) q[118];
swap q[110], q[118];
swap q[65], q[66];
cx q[73], q[66];
cx q[118], q[119];
rz(-pi/4) q[73];
rz(pi/4) q[66];
rz(pi/4) q[118];
rz(-pi/4) q[119];
cx q[73], q[66];
cx q[118], q[119];
cx q[65], q[66];
cx q[118], q[110];
cx q[73], q[66];
cx q[118], q[119];
rz(pi/2) q[118];
sx q[118];
rz(pi/2) q[118];
swap q[110], q[118];
cx q[110], q[100];
rz(-pi/4) q[110];
swap q[100], q[110];
cx q[101], q[100];
rz(pi/4) q[100];
cx q[100], q[110];
rz(-pi/4) q[100];
rz(pi/4) q[110];
cx q[101], q[100];
rz(3*pi/4) q[100];
sx q[100];
rz(pi/2) q[100];
swap q[100], q[101];
cx q[100], q[110];
rz(pi/4) q[100];
rz(-pi/4) q[110];
cx q[100], q[110];
cx q[100], q[101];
cx q[100], q[110];
rz(pi/2) q[100];
sx q[100];
rz(pi/2) q[100];
swap q[100], q[101];
cx q[101], q[102];
rz(-pi/4) q[101];
swap q[101], q[102];
cx q[92], q[102];
rz(pi/4) q[102];
cx q[102], q[101];
rz(-pi/4) q[102];
rz(pi/4) q[101];
cx q[92], q[102];
rz(3*pi/4) q[102];
sx q[102];
rz(pi/2) q[102];
swap q[92], q[102];
cx q[102], q[101];
rz(pi/4) q[102];
rz(-pi/4) q[101];
cx q[102], q[101];
cx q[102], q[92];
cx q[102], q[101];
rz(pi/2) q[102];
sx q[102];
rz(pi/2) q[102];
cx q[102], q[103];
rz(-pi/4) q[102];
swap q[102], q[103];
cx q[103], q[104];
rz(pi/4) q[103];
cx q[103], q[102];
rz(-pi/4) q[103];
rz(pi/4) q[102];
cx q[103], q[104];
rz(3*pi/4) q[103];
sx q[103];
rz(pi/2) q[103];
swap q[103], q[104];
cx q[102], q[103];
rz(-pi/4) q[102];
rz(pi/4) q[103];
cx q[102], q[103];
cx q[104], q[103];
cx q[102], q[103];

// measurement
measure q[77]->c[0];
measure q[94]->c[1];
measure q[96]->c[2];
measure q[98]->c[3];
measure q[63]->c[4];
measure q[58]->c[5];
measure q[53]->c[6];
measure q[40]->c[7];
measure q[43]->c[8];
measure q[47]->c[9];
measure q[64]->c[10];
measure q[69]->c[11];
measure q[48]->c[12];
measure q[68]->c[13];
measure q[74]->c[14];
measure q[89]->c[15];
measure q[84]->c[16];
measure q[72]->c[17];
measure q[81]->c[18];
measure q[78]->c[19];
measure q[65]->c[20];
measure q[85]->c[21];
measure q[87]->c[22];
measure q[106]->c[23];
measure q[104]->c[24];
measure q[92]->c[25];
measure q[100]->c[26];
measure q[118]->c[27];
measure q[75]->c[28];
measure q[90]->c[29];
measure q[95]->c[30];
measure q[97]->c[31];
measure q[60]->c[32];
measure q[59]->c[33];
measure q[41]->c[34];
measure q[39]->c[35];
measure q[45]->c[36];
measure q[46]->c[37];
measure q[54]->c[38];
measure q[67]->c[39];
measure q[50]->c[40];
measure q[55]->c[41];
measure q[70]->c[42];
measure q[88]->c[43];
measure q[82]->c[44];
measure q[62]->c[45];
measure q[80]->c[46];
measure q[79]->c[47];
measure q[73]->c[48];
measure q[86]->c[49];
measure q[93]->c[50];
measure q[105]->c[51];
measure q[102]->c[52];
measure q[101]->c[53];
measure q[110]->c[54];
measure q[119]->c[55];
measure q[76]->c[56];
measure q[61]->c[57];
measure q[44]->c[58];
measure q[49]->c[59];
measure q[83]->c[60];
measure q[66]->c[61];
measure q[103]->c[62];
measure q[117]->c[63];
