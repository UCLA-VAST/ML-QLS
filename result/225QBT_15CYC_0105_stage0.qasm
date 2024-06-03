OPENQASM 2.0;
include "qelib1.inc";
qreg q[225];
creg c[225];
x q[115];
x q[12];
x q[39];
x q[146];
x q[13];
x q[175];
x q[60];
x q[113];
x q[210];
x q[96];
x q[35];
x q[64];
cx q[220], q[219];
cx q[169], q[154];
cx q[56], q[41];
cx q[88], q[103];
cx q[180], q[165];
cx q[7], q[8];
cx q[205], q[204];
cx q[24], q[23];
cx q[95], q[94];
cx q[29], q[44];
cx q[71], q[70];
cx q[208], q[223];
cx q[212], q[197];
cx q[117], q[102];
cx q[19], q[20];
cx q[143], q[142];
cx q[16], q[1];
cx q[173], q[174];
cx q[192], q[207];
cx q[200], q[199];
cx q[85], q[86];
cx q[73], q[74];
cx q[128], q[129];
cx q[49], q[34];
cx q[54], q[53];
cx q[214], q[213];
cx q[4], q[3];
cx q[119], q[134];
cx q[178], q[193];
cx q[125], q[124];
cx q[42], q[57];
cx q[221], q[222];
cx q[203], q[218];
x q[55];
x q[108];
cx q[161], q[176];
x q[21];
x q[150];
x q[155];
x q[219];
x q[205];
x q[73];
x q[54];
x q[57];
x q[146];
cx q[61], q[60];
cx q[125], q[124];
cx q[94], q[93];
cx q[19], q[4];
cx q[115], q[130];
cx q[144], q[143];
cx q[208], q[223];
cx q[12], q[13];
cx q[16], q[1];
cx q[104], q[119];
cx q[85], q[86];
cx q[35], q[36];
cx q[111], q[96];
cx q[56], q[41];
cx q[210], q[211];
cx q[14], q[29];
x q[24];
x q[21];
swap q[163], q[178];
swap q[193], q[194];
swap q[176], q[191];
swap q[103], q[118];
swap q[169], q[184];
swap q[88], q[89];
swap q[217], q[218];
swap q[203], q[204];
swap q[200], q[201];
swap q[173], q[188];
swap q[139], q[154];
swap q[23], q[38];
swap q[206], q[221];
swap q[182], q[197];
swap q[214], q[215];
swap q[2], q[3];
x q[169];
cx q[173], q[158];
cx q[88], q[103];
x q[38];
x q[2];
cx q[182], q[167];
cx q[215], q[214];
cx q[221], q[220];
cx q[139], q[154];
cx q[175], q[176];
cx q[202], q[201];
x q[54];
cx q[15], q[16];
cx q[73], q[74];
cx q[57], q[72];
cx q[19], q[20];
cx q[35], q[36];
x q[29];
x q[210];
swap q[162], q[163];
swap q[194], q[209];
swap q[61], q[62];
swap q[45], q[60];
swap q[130], q[145];
swap q[100], q[115];
swap q[11], q[12];
swap q[190], q[191];
swap q[118], q[133];
swap q[81], q[96];
swap q[111], q[126];
swap q[184], q[185];
swap q[207], q[208];
swap q[86], q[87];
swap q[211], q[212];
swap q[4], q[5];
swap q[0], q[1];
swap q[89], q[104];
x q[12];
cx q[114], q[115];
cx q[76], q[61];
cx q[194], q[193];
cx q[164], q[163];
cx q[111], q[110];
x q[169];
x q[208];
cx q[118], q[117];
cx q[173], q[188];
x q[89];
cx q[220], q[219];
cx q[192], q[207];
cx q[85], q[100];
cx q[0], q[1];
x q[54];
cx q[205], q[190];
cx q[19], q[20];
x q[210];
swap q[147], q[162];
swap q[62], q[63];
swap q[209], q[224];
swap q[145], q[160];
swap q[10], q[11];
swap q[176], q[177];
swap q[126], q[141];
swap q[185], q[186];
swap q[133], q[148];
swap q[202], q[203];
swap q[201], q[216];
swap q[80], q[81];
swap q[57], q[58];
swap q[72], q[87];
swap q[152], q[167];
swap q[213], q[214];
swap q[5], q[6];
swap q[16], q[17];
cx q[5], q[4];
cx q[82], q[81];
cx q[185], q[184];
x q[167];
cx q[114], q[129];
cx q[201], q[200];
cx q[208], q[223];
cx q[74], q[89];
cx q[87], q[86];
cx q[2], q[1];
cx q[12], q[13];
cx q[34], q[19];
x q[210];
swap q[132], q[147];
swap q[9], q[10];
swap q[161], q[162];
swap q[63], q[78];
swap q[186], q[187];
swap q[204], q[205];
swap q[61], q[62];
swap q[75], q[76];
swap q[207], q[222];
swap q[177], q[192];
swap q[84], q[85];
swap q[151], q[152];
swap q[178], q[193];
swap q[149], q[164];
swap q[71], q[72];
swap q[17], q[18];
swap q[58], q[73];
swap q[118], q[133];
cx q[32], q[17];
cx q[4], q[3];
cx q[163], q[178];
cx q[70], q[85];
cx q[183], q[184];
cx q[206], q[205];
cx q[57], q[72];
cx q[162], q[161];
cx q[78], q[93];
cx q[208], q[223];
cx q[24], q[9];
cx q[76], q[75];
swap q[114], q[115];
swap q[128], q[129];
swap q[203], q[204];
swap q[47], q[62];
swap q[81], q[96];
swap q[136], q[151];
swap q[221], q[222];
swap q[87], q[102];
swap q[134], q[149];
swap q[118], q[119];
cx q[166], q[151];
cx q[31], q[32];
cx q[96], q[95];
cx q[16], q[17];
cx q[149], q[164];
cx q[114], q[113];
x q[119];
cx q[118], q[133];
cx q[70], q[85];
cx q[221], q[220];
cx q[5], q[4];
cx q[78], q[93];
x q[24];
swap q[115], q[116];
swap q[160], q[161];
swap q[189], q[204];
swap q[8], q[9];
swap q[162], q[177];
swap q[207], q[208];
swap q[135], q[136];
swap q[183], q[198];
swap q[76], q[91];
swap q[178], q[193];
swap q[148], q[163];
swap q[42], q[57];
cx q[166], q[165];
x q[164];
x q[116];
cx q[199], q[198];
cx q[189], q[188];
cx q[128], q[113];
cx q[183], q[182];
x q[160];
x q[177];
cx q[161], q[146];
cx q[17], q[18];
cx q[95], q[94];
cx q[194], q[193];
cx q[9], q[10];
cx q[86], q[85];
cx q[118], q[133];
swap q[204], q[205];
swap q[206], q[207];
swap q[96], q[97];
swap q[63], q[78];
swap q[27], q[42];
swap q[221], q[222];
swap q[5], q[6];
swap q[30], q[31];
swap q[56], q[57];
swap q[70], q[71];
swap q[136], q[151];
swap q[134], q[149];
cx q[43], q[42];
cx q[46], q[31];
cx q[137], q[136];
x q[5];
cx q[166], q[165];
cx q[134], q[149];
cx q[183], q[182];
cx q[173], q[188];
cx q[72], q[71];
cx q[151], q[150];
x q[116];
x q[146];
cx q[129], q[128];
cx q[94], q[93];
cx q[221], q[220];
x q[9];
x q[133];
swap q[203], q[204];
swap q[205], q[206];
swap q[197], q[198];
swap q[207], q[208];
swap q[199], q[200];
swap q[48], q[63];
swap q[80], q[95];
swap q[179], q[194];
swap q[222], q[223];
swap q[15], q[30];
cx q[206], q[207];
cx q[213], q[198];
x q[137];
cx q[167], q[166];
cx q[119], q[134];
cx q[15], q[16];
x q[95];
cx q[57], q[72];
cx q[5], q[6];
x q[9];
swap q[113], q[128];
swap q[129], q[130];
swap q[173], q[174];
swap q[204], q[205];
swap q[202], q[203];
swap q[208], q[209];
swap q[33], q[48];
swap q[94], q[109];
swap q[92], q[93];
swap q[219], q[220];
swap q[27], q[42];
swap q[43], q[58];
x q[94];
x q[93];
cx q[158], q[173];
cx q[203], q[218];
x q[27];
cx q[202], q[217];
x q[220];
cx q[198], q[197];
cx q[166], q[165];
cx q[213], q[212];
cx q[16], q[17];
cx q[137], q[136];
cx q[134], q[149];
cx q[92], q[91];
cx q[9], q[24];
swap q[130], q[145];
swap q[98], q[113];
swap q[189], q[204];
swap q[109], q[110];
swap q[193], q[208];
swap q[4], q[5];
swap q[6], q[7];
swap q[56], q[57];
swap q[33], q[34];
swap q[72], q[73];
swap q[104], q[119];
x q[93];
x q[94];
x q[6];
cx q[58], q[57];
x q[203];
cx q[202], q[217];
cx q[27], q[42];
x q[119];
cx q[208], q[223];
cx q[109], q[108];
x q[145];
x q[165];
cx q[136], q[151];
cx q[198], q[183];
cx q[219], q[220];
x q[149];
cx q[111], q[110];
swap q[24], q[25];
swap q[8], q[9];
swap q[90], q[91];
swap q[77], q[92];
swap q[3], q[4];
swap q[15], q[16];
swap q[41], q[56];
cx q[106], q[91];
cx q[57], q[72];
cx q[94], q[93];
x q[56];
cx q[27], q[42];
cx q[104], q[119];
cx q[202], q[201];
cx q[15], q[30];
cx q[21], q[6];
cx q[59], q[58];
cx q[160], q[145];
cx q[109], q[108];
cx q[213], q[198];
cx q[219], q[220];
swap q[216], q[217];
swap q[62], q[77];
swap q[193], q[208];
swap q[90], q[105];
swap q[3], q[18];
swap q[223], q[224];
swap q[26], q[41];
swap q[150], q[151];
cx q[40], q[41];
cx q[218], q[217];
x q[223];
cx q[17], q[18];
cx q[165], q[150];
cx q[27], q[26];
x q[109];
cx q[213], q[212];
cx q[219], q[220];
swap q[130], q[145];
swap q[160], q[175];
swap q[200], q[201];
swap q[62], q[63];
swap q[107], q[108];
swap q[76], q[91];
swap q[2], q[3];
swap q[72], q[73];
swap q[42], q[57];
swap q[6], q[7];
swap q[44], q[59];
cx q[40], q[41];
cx q[77], q[62];
x q[218];
x q[76];
cx q[57], q[72];
cx q[6], q[5];
cx q[2], q[3];
x q[165];
x q[18];
cx q[150], q[151];
cx q[200], q[199];
x q[219];
swap q[115], q[130];
swap q[159], q[160];
swap q[175], q[190];
swap q[197], q[212];
swap q[213], q[214];
swap q[202], q[217];
swap q[48], q[63];
swap q[106], q[107];
swap q[16], q[17];
swap q[11], q[26];
swap q[27], q[42];
swap q[73], q[88];
x q[130];
cx q[55], q[40];
x q[218];
cx q[175], q[174];
x q[202];
cx q[17], q[32];
cx q[15], q[16];
cx q[57], q[72];
x q[6];
x q[165];
cx q[150], q[151];
cx q[200], q[199];
swap q[144], q[159];
swap q[196], q[197];
swap q[106], q[121];
swap q[1], q[2];
swap q[3], q[4];
swap q[26], q[27];
swap q[73], q[74];
x q[197];
cx q[28], q[27];
x q[130];
cx q[59], q[74];
cx q[158], q[159];
cx q[2], q[3];
cx q[73], q[88];
cx q[55], q[40];
cx q[203], q[218];
cx q[4], q[5];
x q[6];
cx q[42], q[57];
cx q[165], q[150];
cx q[17], q[32];
swap q[175], q[176];
swap q[173], q[174];
swap q[195], q[196];
swap q[199], q[214];
swap q[120], q[121];
swap q[15], q[30];
swap q[11], q[26];
x q[175];
cx q[196], q[181];
x q[27];
x q[158];
x q[173];
cx q[177], q[176];
cx q[55], q[40];
x q[5];
x q[165];
swap q[203], q[204];
swap q[217], q[218];
swap q[213], q[214];
swap q[30], q[45];
swap q[11], q[12];
swap q[41], q[42];
swap q[88], q[103];
swap q[73], q[74];
swap q[3], q[4];
cx q[160], q[175];
cx q[28], q[27];
cx q[197], q[196];
cx q[215], q[214];
cx q[73], q[88];
x q[158];
cx q[10], q[11];
cx q[59], q[74];
x q[218];
x q[4];
cx q[18], q[3];
swap q[212], q[213];
swap q[176], q[191];
swap q[45], q[60];
swap q[103], q[118];
x q[176];
x q[215];
cx q[28], q[27];
cx q[88], q[103];
cx q[173], q[158];
cx q[10], q[9];
x q[11];
x q[3];
x q[4];
x q[212];
swap q[145], q[160];
swap q[199], q[214];
swap q[60], q[75];
swap q[181], q[196];
swap q[73], q[74];
cx q[211], q[196];
cx q[145], q[130];
cx q[176], q[161];
cx q[27], q[42];
cx q[216], q[215];
cx q[74], q[89];
x q[10];
x q[214];
x q[4];
swap q[173], q[188];
swap q[143], q[158];
swap q[75], q[90];
swap q[103], q[118];
swap q[8], q[9];
cx q[181], q[196];
cx q[216], q[215];
cx q[74], q[89];
cx q[161], q[146];
cx q[10], q[11];
cx q[118], q[133];
x q[9];
cx q[5], q[4];
swap q[142], q[143];
swap q[145], q[160];
swap q[115], q[130];
swap q[188], q[189];
swap q[175], q[176];
swap q[12], q[27];
cx q[74], q[89];
cx q[13], q[12];
cx q[27], q[26];
cx q[190], q[189];
cx q[118], q[133];
cx q[215], q[214];
x q[4];
cx q[6], q[5];
swap q[141], q[142];
swap q[114], q[115];
swap q[173], q[188];
swap q[216], q[217];
swap q[196], q[197];
cx q[196], q[181];
cx q[116], q[115];
x q[216];
cx q[218], q[217];
cx q[13], q[12];
cx q[200], q[215];
swap q[126], q[141];
swap q[99], q[114];
swap q[158], q[173];
swap q[190], q[191];
swap q[103], q[118];
swap q[88], q[89];
cx q[196], q[181];
cx q[89], q[104];
cx q[73], q[88];
x q[190];
cx q[217], q[216];
x q[13];
x q[218];
cx q[215], q[214];
swap q[98], q[99];
swap q[143], q[158];
swap q[114], q[115];
x q[104];
cx q[74], q[89];
x q[88];
cx q[218], q[217];
cx q[13], q[12];
cx q[190], q[189];
swap q[83], q[98];
swap q[142], q[143];
swap q[200], q[215];
swap q[213], q[214];
swap q[196], q[211];
swap q[181], q[182];
x q[196];
cx q[128], q[143];
cx q[89], q[104];
cx q[210], q[211];
cx q[88], q[103];
cx q[217], q[216];
cx q[13], q[12];
cx q[212], q[213];
swap q[82], q[83];
swap q[127], q[142];
swap q[200], q[201];
swap q[190], q[191];
swap q[182], q[183];
swap q[59], q[74];
cx q[157], q[142];
cx q[199], q[200];
cx q[128], q[143];
cx q[196], q[197];
cx q[74], q[59];
cx q[210], q[211];
cx q[88], q[103];
swap q[81], q[82];
swap q[112], q[127];
swap q[83], q[84];
swap q[217], q[218];
swap q[215], q[216];
swap q[191], q[192];
swap q[13], q[28];
cx q[142], q[141];
cx q[191], q[206];
cx q[199], q[200];
cx q[13], q[14];
cx q[112], q[97];
cx q[128], q[143];
x q[197];
cx q[210], q[195];
x q[216];
x q[218];
swap q[66], q[81];
swap q[82], q[83];
swap q[157], q[172];
swap q[28], q[43];
swap q[103], q[118];
cx q[173], q[172];
cx q[158], q[157];
x q[28];
cx q[13], q[14];
x q[199];
cx q[196], q[197];
cx q[217], q[218];
cx q[216], q[201];
swap q[65], q[66];
swap q[81], q[82];
swap q[112], q[113];
swap q[206], q[207];
swap q[43], q[58];
cx q[112], q[127];
cx q[51], q[66];
x q[28];
cx q[207], q[206];
cx q[42], q[43];
cx q[199], q[198];
cx q[14], q[29];
swap q[113], q[114];
swap q[157], q[172];
swap q[12], q[13];
cx q[66], q[65];
cx q[157], q[172];
cx q[113], q[128];
x q[12];
cx q[44], q[43];
cx q[28], q[27];
swap q[51], q[52];
swap q[114], q[115];
swap q[205], q[206];
swap q[126], q[127];
swap q[41], q[42];
swap q[199], q[214];
x q[126];
cx q[66], q[81];
cx q[200], q[199];
cx q[207], q[206];
cx q[128], q[143];
cx q[40], q[41];
x q[214];
cx q[11], q[12];
cx q[29], q[28];
swap q[50], q[65];
swap q[115], q[130];
swap q[52], q[53];
swap q[171], q[172];
swap q[42], q[43];
cx q[65], q[64];
cx q[141], q[126];
x q[53];
x q[50];
x q[199];
x q[200];
x q[214];
cx q[29], q[28];
cx q[10], q[11];
cx q[42], q[27];
swap q[81], q[96];
swap q[51], q[66];
swap q[130], q[131];
swap q[205], q[206];
swap q[207], q[222];
swap q[170], q[171];
swap q[127], q[128];
swap q[26], q[41];
cx q[171], q[172];
cx q[66], q[81];
cx q[130], q[129];
cx q[126], q[125];
x q[207];
cx q[52], q[53];
x q[64];
cx q[51], q[50];
cx q[223], q[222];
x q[199];
cx q[205], q[204];
x q[28];
cx q[14], q[29];
cx q[9], q[10];
cx q[11], q[12];
x q[42];
cx q[215], q[214];
swap q[96], q[97];
swap q[131], q[132];
swap q[141], q[156];
cx q[129], q[144];
x q[207];
cx q[38], q[53];
cx q[198], q[199];
x q[11];
cx q[42], q[43];
swap q[97], q[98];
swap q[80], q[81];
swap q[37], q[52];
swap q[115], q[130];
swap q[223], q[224];
swap q[221], q[222];
swap q[171], q[186];
swap q[8], q[9];
swap q[12], q[13];
swap q[29], q[44];
x q[171];
cx q[186], q[185];
cx q[66], q[81];
x q[115];
x q[207];
x q[29];
cx q[200], q[199];
cx q[13], q[14];
cx q[10], q[11];
swap q[98], q[99];
swap q[52], q[53];
swap q[23], q[38];
swap q[129], q[130];
swap q[222], q[223];
cx q[185], q[184];
cx q[130], q[131];
x q[222];
cx q[96], q[81];
cx q[52], q[51];
x q[29];
swap q[23], q[24];
swap q[98], q[113];
swap q[9], q[10];
swap q[11], q[12];
swap q[200], q[201];
cx q[22], q[23];
cx q[186], q[185];
cx q[130], q[131];
x q[222];
cx q[96], q[111];
cx q[14], q[29];
cx q[12], q[13];
swap q[81], q[82];
swap q[50], q[51];
swap q[83], q[98];
cx q[98], q[97];
cx q[23], q[38];
cx q[185], q[184];
x q[131];
cx q[115], q[130];
x q[222];
cx q[50], q[51];
cx q[14], q[29];
swap q[81], q[96];
swap q[110], q[111];
cx q[97], q[112];
cx q[22], q[23];
cx q[39], q[38];
cx q[130], q[131];
cx q[222], q[207];
cx q[110], q[111];
swap q[81], q[82];
swap q[51], q[52];
swap q[114], q[115];
swap q[49], q[50];
swap q[183], q[184];
x q[115];
cx q[98], q[97];
cx q[65], q[50];
cx q[23], q[38];
x q[207];
cx q[110], q[109];
swap q[52], q[53];
swap q[131], q[132];
swap q[112], q[127];
cx q[142], q[127];
cx q[100], q[115];
cx q[65], q[50];
x q[38];
cx q[98], q[97];
cx q[222], q[207];
swap q[108], q[109];
swap q[132], q[147];
swap q[53], q[68];
swap q[131], q[146];
cx q[117], q[132];
cx q[124], q[109];
cx q[100], q[115];
cx q[50], q[35];
cx q[39], q[38];
cx q[131], q[116];
x q[222];
swap q[107], q[108];
swap q[147], q[162];
swap q[68], q[69];
swap q[142], q[157];
swap q[126], q[127];
swap q[206], q[207];
x q[126];
x q[132];
cx q[124], q[109];
x q[142];
x q[157];
x q[147];
cx q[23], q[38];
cx q[39], q[54];
cx q[100], q[115];
swap q[106], q[107];
swap q[162], q[177];
swap q[69], q[84];
swap q[131], q[146];
swap q[101], q[116];
swap q[49], q[50];
cx q[68], q[69];
cx q[107], q[92];
cx q[108], q[109];
cx q[117], q[132];
cx q[126], q[127];
x q[124];
x q[147];
cx q[142], q[157];
x q[50];
cx q[100], q[115];
swap q[177], q[178];
swap q[91], q[106];
swap q[84], q[99];
swap q[146], q[161];
swap q[86], q[101];
swap q[48], q[49];
swap q[24], q[39];
swap q[22], q[23];
cx q[68], q[67];
cx q[93], q[92];
cx q[126], q[127];
cx q[107], q[108];
x q[23];
cx q[117], q[132];
cx q[125], q[124];
cx q[7], q[22];
cx q[51], q[50];
swap q[178], q[193];
swap q[84], q[85];
swap q[142], q[143];
swap q[86], q[87];
swap q[101], q[116];
swap q[99], q[114];
cx q[86], q[101];
cx q[94], q[93];
cx q[107], q[108];
x q[23];
x q[143];
cx q[116], q[117];
cx q[132], q[147];
x q[22];
swap q[193], q[208];
swap q[91], q[92];
swap q[68], q[69];
swap q[66], q[67];
swap q[111], q[126];
swap q[163], q[178];
swap q[72], q[87];
x q[178];
cx q[86], q[101];
cx q[69], q[68];
cx q[91], q[76];
cx q[66], q[67];
x q[72];
cx q[117], q[132];
cx q[23], q[24];
cx q[142], q[143];
cx q[7], q[22];
swap q[208], q[209];
swap q[148], q[163];
cx q[178], q[193];
x q[208];
cx q[86], q[101];
cx q[164], q[163];
cx q[69], q[68];
cx q[209], q[224];
cx q[23], q[24];
cx q[117], q[132];
cx q[57], q[72];
cx q[7], q[8];
swap q[61], q[76];
swap q[66], q[81];
swap q[21], q[22];
cx q[75], q[76];
x q[193];
cx q[164], q[163];
cx q[68], q[67];
cx q[208], q[207];
x q[22];
cx q[209], q[224];
cx q[8], q[9];
swap q[46], q[61];
swap q[81], q[96];
swap q[69], q[84];
swap q[72], q[87];
swap q[101], q[116];
swap q[131], q[132];
cx q[164], q[179];
cx q[178], q[193];
x q[46];
x q[132];
cx q[67], q[66];
x q[22];
x q[101];
x q[72];
swap q[61], q[62];
swap q[84], q[85];
swap q[68], q[83];
swap q[96], q[97];
swap q[87], q[102];
swap q[54], q[69];
swap q[148], q[163];
swap q[194], q[209];
x q[209];
x q[62];
cx q[148], q[163];
cx q[95], q[96];
cx q[86], q[85];
cx q[178], q[193];
cx q[83], q[82];
cx q[22], q[23];
swap q[66], q[81];
swap q[52], q[67];
swap q[53], q[68];
swap q[102], q[117];
swap q[54], q[55];
swap q[179], q[194];
cx q[67], q[66];
cx q[61], q[62];
cx q[96], q[111];
x q[209];
cx q[39], q[54];
cx q[149], q[148];
cx q[101], q[102];
cx q[23], q[24];
swap q[37], q[52];
swap q[83], q[98];
swap q[84], q[85];
swap q[116], q[117];
swap q[163], q[178];
cx q[67], q[66];
cx q[95], q[96];
cx q[85], q[70];
cx q[148], q[147];
cx q[117], q[132];
cx q[23], q[24];
swap q[61], q[76];
swap q[111], q[126];
swap q[87], q[102];
swap q[100], q[101];
swap q[39], q[40];
swap q[134], q[149];
cx q[67], q[66];
cx q[61], q[76];
cx q[103], q[102];
cx q[149], q[164];
x q[126];
cx q[85], q[70];
cx q[147], q[132];
cx q[54], q[39];
cx q[22], q[23];
swap q[96], q[111];
swap q[86], q[101];
swap q[148], q[163];
swap q[87], q[88];
x q[67];
x q[66];
cx q[103], q[102];
cx q[149], q[164];
cx q[111], q[126];
cx q[81], q[96];
cx q[147], q[132];
cx q[22], q[7];
swap q[76], q[77];
swap q[39], q[40];
swap q[53], q[54];
swap q[85], q[86];
swap q[55], q[70];
cx q[71], q[70];
cx q[67], q[66];
cx q[76], q[77];
cx q[103], q[102];
cx q[82], q[81];
cx q[86], q[85];
swap q[126], q[127];
swap q[95], q[96];
swap q[52], q[53];
swap q[25], q[40];
swap q[131], q[132];
swap q[134], q[149];
x q[53];
cx q[71], q[70];
cx q[134], q[149];
cx q[40], q[39];
cx q[132], q[117];
cx q[110], q[95];
x q[85];
x q[86];
cx q[25], q[26];
cx q[87], q[102];
swap q[75], q[76];
swap q[127], q[142];
swap q[82], q[83];
swap q[80], q[81];
swap q[51], q[52];
swap q[116], q[131];
swap q[88], q[103];
x q[76];
cx q[56], q[71];
cx q[38], q[53];
x q[142];
cx q[73], q[88];
x q[132];
cx q[51], q[50];
x q[26];
x q[85];
x q[86];
swap q[126], q[127];
swap q[94], q[95];
swap q[110], q[111];
swap q[83], q[84];
swap q[79], q[80];
swap q[148], q[149];
swap q[119], q[134];
swap q[101], q[102];
cx q[80], q[81];
cx q[76], q[61];
x q[102];
x q[149];
cx q[41], q[56];
x q[127];
cx q[104], q[119];
x q[88];
cx q[38], q[39];
cx q[133], q[148];
cx q[74], q[73];
cx q[11], q[26];
cx q[132], q[117];
x q[86];
swap q[111], q[126];
swap q[84], q[99];
swap q[68], q[83];
swap q[51], q[52];
x q[83];
cx q[91], q[76];
cx q[69], q[84];
cx q[134], q[149];
cx q[41], q[56];
x q[102];
x q[51];
x q[88];
cx q[40], q[39];
cx q[142], q[127];
cx q[118], q[133];
x q[132];
cx q[25], q[26];
swap q[126], q[141];
swap q[80], q[95];
swap q[66], q[81];
swap q[67], q[68];
swap q[59], q[74];
x q[95];
x q[66];
cx q[83], q[82];
x q[81];
cx q[76], q[75];
cx q[74], q[89];
cx q[103], q[102];
cx q[51], q[50];
cx q[40], q[39];
cx q[73], q[88];
swap q[91], q[106];
swap q[141], q[156];
swap q[125], q[126];
swap q[68], q[69];
swap q[149], q[164];
swap q[55], q[56];
swap q[119], q[134];
swap q[26], q[41];
cx q[140], q[125];
x q[95];
x q[83];
cx q[107], q[106];
cx q[81], q[96];
cx q[75], q[90];
cx q[74], q[89];
cx q[27], q[26];
cx q[103], q[102];
cx q[51], q[50];
cx q[73], q[88];
swap q[91], q[92];
swap q[156], q[157];
swap q[53], q[68];
swap q[164], q[179];
swap q[69], q[70];
swap q[38], q[39];
cx q[140], q[125];
cx q[95], q[80];
x q[70];
cx q[83], q[82];
cx q[107], q[108];
cx q[53], q[68];
x q[39];
cx q[89], q[104];
x q[74];
cx q[27], q[26];
cx q[102], q[117];
swap q[90], q[105];
swap q[60], q[75];
swap q[92], q[93];
swap q[91], q[106];
swap q[157], q[158];
cx q[155], q[140];
cx q[95], q[80];
cx q[91], q[76];
cx q[71], q[70];
x q[106];
cx q[74], q[89];
cx q[39], q[38];
swap q[45], q[60];
swap q[90], q[105];
swap q[158], q[159];
swap q[52], q[53];
swap q[81], q[82];
swap q[68], q[83];
swap q[104], q[119];
cx q[158], q[173];
cx q[141], q[140];
x q[105];
x q[95];
cx q[82], q[81];
cx q[46], q[45];
cx q[91], q[76];
cx q[72], q[71];
cx q[119], q[104];
x q[106];
cx q[59], q[74];
swap q[159], q[160];
swap q[79], q[80];
swap q[37], q[52];
swap q[155], q[170];
swap q[83], q[98];
swap q[24], q[39];
x q[155];
x q[170];
cx q[110], q[95];
x q[76];
cx q[57], q[72];
cx q[119], q[134];
swap q[46], q[47];
swap q[30], q[45];
swap q[160], q[175];
swap q[173], q[174];
swap q[143], q[158];
swap q[126], q[141];
swap q[82], q[83];
swap q[37], q[38];
cx q[173], q[172];
cx q[45], q[60];
x q[155];
cx q[61], q[46];
x q[82];
x q[170];
cx q[140], q[141];
cx q[128], q[143];
cx q[57], q[72];
cx q[39], q[38];
cx q[134], q[149];
swap q[32], q[47];
swap q[175], q[176];
swap q[80], q[95];
swap q[160], q[161];
swap q[110], q[111];
x q[175];
cx q[173], q[172];
cx q[45], q[60];
cx q[32], q[17];
swap q[176], q[191];
swap q[79], q[80];
swap q[159], q[160];
swap q[126], q[141];
swap q[96], q[111];
x q[175];
cx q[145], q[160];
cx q[176], q[177];
cx q[173], q[174];
cx q[156], q[141];
cx q[60], q[75];
x q[80];
cx q[206], q[191];
cx q[81], q[96];
swap q[32], q[47];
swap q[30], q[45];
swap q[172], q[187];
cx q[171], q[172];
cx q[145], q[144];
cx q[160], q[175];
x q[30];
cx q[176], q[177];
cx q[173], q[174];
cx q[60], q[61];
x q[32];
swap q[47], q[62];
swap q[45], q[46];
swap q[206], q[221];
swap q[190], q[191];
swap q[186], q[187];
cx q[157], q[172];
cx q[187], q[186];
cx q[177], q[162];
cx q[173], q[174];
cx q[17], q[32];
cx q[221], q[220];
swap q[62], q[77];
swap q[45], q[60];
swap q[130], q[145];
swap q[176], q[191];
cx q[171], q[172];
x q[177];
cx q[75], q[60];
cx q[61], q[62];
cx q[30], q[45];
cx q[219], q[220];
swap q[77], q[92];
swap q[145], q[160];
swap q[156], q[157];
cx q[130], q[145];
x q[156];
cx q[186], q[171];
x q[77];
swap q[91], q[92];
swap q[45], q[46];
swap q[75], q[90];
x q[145];
cx q[156], q[155];
x q[186];
cx q[105], q[90];
cx q[46], q[61];
cx q[77], q[92];
x q[45];
swap q[129], q[130];
cx q[129], q[144];
cx q[156], q[155];
cx q[120], q[105];
cx q[76], q[61];
cx q[92], q[91];
swap q[77], q[78];
swap q[46], q[47];
cx q[129], q[144];
cx q[156], q[155];
cx q[105], q[90];
cx q[107], q[92];
swap q[62], q[77];
cx q[62], q[63];
x q[156];
x q[105];
swap q[107], q[108];
swap q[91], q[92];
swap q[144], q[159];
swap q[114], q[129];
swap q[108], q[123];
swap q[76], q[91];
swap q[159], q[160];
cx q[76], q[91];
swap q[123], q[138];
swap q[160], q[175];
swap q[158], q[159];
cx q[160], q[161];
cx q[174], q[159];
swap q[138], q[153];
swap q[61], q[76];
swap q[175], q[190];
cx q[138], q[123];
x q[160];
cx q[161], q[146];
cx q[190], q[191];
cx q[173], q[174];
cx q[76], q[77];
swap q[46], q[61];
swap q[158], q[159];
cx q[145], q[160];
cx q[191], q[206];
cx q[158], q[157];
x q[76];
swap q[31], q[46];
swap q[108], q[123];
swap q[173], q[188];
cx q[123], q[122];
cx q[108], q[93];
x q[46];
cx q[160], q[161];
cx q[206], q[207];
swap q[16], q[31];
swap q[187], q[188];
swap q[130], q[145];
swap q[191], q[192];
cx q[61], q[46];
x q[145];
cx q[161], q[162];
cx q[187], q[186];
cx q[208], q[207];
cx q[192], q[177];
x q[191];
swap q[1], q[16];
swap q[78], q[93];
swap q[108], q[109];
swap q[123], q[124];
swap q[160], q[175];
swap q[115], q[130];
x q[109];
cx q[78], q[63];
cx q[139], q[124];
cx q[61], q[46];
cx q[123], q[108];
cx q[190], q[175];
cx q[144], q[145];
x q[177];
swap q[0], q[1];
swap q[187], q[188];
swap q[160], q[161];
swap q[129], q[130];
swap q[185], q[186];
swap q[147], q[162];
swap q[206], q[207];
swap q[208], q[223];
x q[1];
x q[186];
x q[61];
cx q[115], q[130];
cx q[129], q[144];
cx q[192], q[207];
cx q[224], q[223];
swap q[31], q[46];
swap q[48], q[63];
swap q[139], q[154];
swap q[123], q[124];
swap q[188], q[189];
swap q[159], q[160];
swap q[193], q[208];
swap q[205], q[206];
cx q[139], q[138];
cx q[123], q[122];
x q[1];
cx q[169], q[154];
cx q[78], q[63];
cx q[178], q[193];
cx q[186], q[171];
cx q[208], q[209];
x q[61];
x q[129];
cx q[192], q[207];
cx q[224], q[223];
swap q[33], q[48];
swap q[173], q[188];
swap q[160], q[175];
swap q[124], q[125];
x q[124];
cx q[123], q[122];
x q[139];
cx q[1], q[2];
cx q[203], q[188];
x q[169];
cx q[78], q[63];
cx q[194], q[209];
cx q[163], q[178];
cx q[60], q[61];
x q[193];
x q[175];
cx q[207], q[206];
swap q[160], q[161];
swap q[171], q[172];
cx q[124], q[109];
x q[171];
cx q[188], q[187];
x q[2];
x q[169];
cx q[16], q[1];
cx q[203], q[202];
cx q[63], q[48];
x q[163];
x q[160];
cx q[175], q[190];
swap q[60], q[75];
swap q[46], q[61];
swap q[78], q[93];
swap q[146], q[161];
swap q[121], q[122];
swap q[123], q[138];
swap q[208], q[209];
swap q[179], q[194];
x q[122];
cx q[138], q[123];
x q[171];
x q[161];
cx q[124], q[125];
cx q[202], q[187];
x q[179];
cx q[194], q[209];
cx q[3], q[2];
cx q[163], q[162];
cx q[204], q[203];
x q[75];
x q[160];
cx q[208], q[193];
cx q[175], q[190];
swap q[46], q[47];
swap q[48], q[49];
swap q[93], q[108];
swap q[15], q[16];
cx q[121], q[122];
x q[108];
cx q[16], q[31];
cx q[64], q[49];
cx q[172], q[171];
x q[202];
x q[15];
cx q[161], q[176];
cx q[188], q[187];
cx q[126], q[125];
cx q[194], q[179];
cx q[163], q[178];
cx q[205], q[204];
cx q[90], q[75];
swap q[47], q[62];
swap q[138], q[153];
swap q[147], q[162];
cx q[154], q[153];
x q[171];
cx q[16], q[15];
cx q[202], q[187];
cx q[161], q[176];
x q[162];
x q[62];
cx q[173], q[188];
cx q[194], q[209];
cx q[204], q[203];
cx q[148], q[163];
x q[178];
cx q[105], q[90];
swap q[63], q[64];
swap q[34], q[49];
swap q[122], q[123];
x q[49];
x q[122];
cx q[154], q[139];
cx q[171], q[170];
x q[15];
cx q[202], q[201];
x q[162];
x q[173];
cx q[77], q[62];
cx q[187], q[186];
cx q[188], q[189];
x q[204];
swap q[63], q[78];
swap q[64], q[79];
swap q[16], q[31];
swap q[34], q[35];
swap q[146], q[161];
swap q[152], q[153];
swap q[133], q[148];
swap q[105], q[120];
cx q[153], q[168];
cx q[122], q[123];
cx q[49], q[48];
cx q[16], q[1];
cx q[154], q[139];
cx q[93], q[78];
x q[171];
x q[105];
cx q[202], q[201];
cx q[177], q[162];
cx q[188], q[187];
cx q[148], q[133];
swap q[35], q[36];
swap q[145], q[146];
swap q[137], q[152];
cx q[153], q[168];
x q[152];
cx q[35], q[20];
cx q[123], q[108];
x q[154];
x q[16];
cx q[1], q[2];
cx q[145], q[144];
cx q[202], q[201];
cx q[156], q[171];
cx q[192], q[177];
cx q[106], q[105];
swap q[47], q[48];
swap q[49], q[64];
swap q[146], q[161];
swap q[63], q[78];
swap q[148], q[163];
cx q[153], q[168];
cx q[48], q[33];
x q[35];
cx q[19], q[20];
cx q[152], q[137];
x q[47];
cx q[169], q[154];
cx q[31], q[16];
cx q[176], q[161];
cx q[145], q[144];
cx q[201], q[186];
swap q[64], q[65];
swap q[122], q[123];
swap q[93], q[108];
swap q[163], q[164];
swap q[2], q[3];
swap q[171], q[172];
x q[64];
cx q[33], q[34];
cx q[19], q[20];
cx q[152], q[137];
x q[171];
cx q[66], q[65];
x q[2];
cx q[107], q[108];
x q[163];
cx q[145], q[144];
swap q[30], q[31];
swap q[167], q[168];
swap q[93], q[94];
swap q[139], q[154];
swap q[175], q[176];
swap q[185], q[186];
cx q[183], q[168];
cx q[48], q[33];
cx q[34], q[35];
x q[154];
cx q[64], q[79];
cx q[109], q[94];
x q[171];
cx q[31], q[46];
cx q[30], q[15];
cx q[2], q[3];
cx q[107], q[92];
x q[175];
x q[186];
cx q[129], q[144];
swap q[145], q[146];
swap q[136], q[137];
swap q[166], q[167];
swap q[124], q[139];
x q[167];
x q[137];
cx q[153], q[168];
cx q[33], q[34];
cx q[94], q[79];
cx q[140], q[139];
cx q[47], q[46];
cx q[3], q[4];
x q[171];
cx q[107], q[108];
swap q[0], q[15];
swap q[144], q[159];
swap q[129], q[130];
swap q[121], q[136];
swap q[35], q[36];
swap q[124], q[125];
x q[167];
x q[136];
cx q[48], q[33];
cx q[183], q[168];
cx q[129], q[144];
cx q[154], q[139];
x q[0];
cx q[31], q[46];
x q[3];
cx q[5], q[4];
swap q[47], q[62];
swap q[159], q[174];
swap q[152], q[153];
swap q[21], q[36];
swap q[130], q[131];
swap q[109], q[124];
swap q[64], q[79];
swap q[34], q[35];
swap q[140], q[141];
cx q[137], q[136];
cx q[155], q[140];
cx q[123], q[124];
cx q[159], q[158];
cx q[141], q[126];
cx q[129], q[144];
cx q[48], q[33];
x q[0];
x q[47];
swap q[30], q[31];
swap q[46], q[61];
swap q[183], q[184];
swap q[153], q[168];
swap q[130], q[145];
swap q[62], q[63];
swap q[94], q[109];
swap q[20], q[21];
swap q[19], q[34];
swap q[36], q[51];
cx q[121], q[136];
cx q[115], q[130];
cx q[153], q[168];
cx q[155], q[170];
cx q[94], q[79];
cx q[199], q[184];
cx q[35], q[20];
x q[51];
cx q[123], q[124];
cx q[127], q[126];
x q[144];
cx q[61], q[76];
cx q[140], q[141];
swap q[31], q[46];
swap q[114], q[129];
swap q[48], q[63];
swap q[36], q[37];
swap q[157], q[158];
swap q[18], q[33];
x q[31];
cx q[137], q[136];
cx q[169], q[168];
x q[35];
cx q[20], q[21];
cx q[155], q[170];
x q[51];
cx q[36], q[37];
swap q[60], q[61];
swap q[76], q[77];
swap q[99], q[114];
swap q[78], q[79];
swap q[127], q[142];
swap q[111], q[126];
swap q[130], q[145];
swap q[122], q[123];
swap q[124], q[139];
x q[31];
cx q[136], q[151];
x q[169];
x q[35];
x q[21];
cx q[124], q[139];
cx q[155], q[170];
cx q[145], q[144];
swap q[46], q[61];
swap q[114], q[129];
swap q[99], q[100];
swap q[63], q[78];
swap q[110], q[111];
swap q[142], q[157];
swap q[60], q[75];
cx q[112], q[111];
x q[31];
cx q[121], q[136];
cx q[45], q[46];
cx q[150], q[151];
cx q[63], q[62];
x q[157];
x q[61];
cx q[184], q[169];
cx q[115], q[100];
x q[60];
cx q[21], q[6];
swap q[84], q[99];
swap q[78], q[93];
swap q[109], q[124];
swap q[129], q[130];
cx q[99], q[114];
cx q[112], q[111];
cx q[123], q[124];
x q[129];
cx q[63], q[64];
cx q[46], q[31];
cx q[121], q[122];
cx q[165], q[150];
x q[157];
cx q[100], q[85];
cx q[76], q[61];
swap q[94], q[109];
swap q[69], q[84];
swap q[93], q[108];
cx q[99], q[114];
cx q[113], q[112];
cx q[111], q[126];
x q[109];
cx q[123], q[108];
swap q[165], q[180];
swap q[135], q[150];
swap q[100], q[115];
swap q[49], q[64];
swap q[63], q[78];
swap q[83], q[84];
swap q[69], q[70];
swap q[45], q[46];
x q[165];
x q[150];
cx q[83], q[68];
x q[114];
x q[69];
x q[63];
cx q[112], q[111];
x q[46];
cx q[110], q[109];
x q[100];
x q[45];
x q[108];
cx q[123], q[124];
swap q[115], q[116];
swap q[55], q[70];
swap q[34], q[49];
swap q[113], q[128];
swap q[135], q[136];
swap q[98], q[99];
cx q[166], q[165];
x q[70];
x q[98];
cx q[68], q[67];
x q[55];
cx q[69], q[84];
cx q[115], q[130];
cx q[99], q[114];
cx q[128], q[143];
cx q[63], q[48];
x q[46];
cx q[34], q[35];
cx q[94], q[109];
swap q[111], q[126];
x q[165];
cx q[181], q[166];
cx q[83], q[68];
cx q[70], q[71];
cx q[67], q[82];
cx q[99], q[114];
cx q[129], q[128];
cx q[115], q[130];
cx q[111], q[96];
cx q[63], q[48];
swap q[54], q[69];
swap q[143], q[158];
swap q[34], q[49];
x q[166];
cx q[150], q[165];
cx q[182], q[181];
x q[69];
cx q[19], q[34];
cx q[70], q[55];
cx q[142], q[143];
cx q[158], q[159];
cx q[68], q[67];
swap q[82], q[97];
swap q[130], q[145];
x q[182];
cx q[181], q[166];
cx q[150], q[165];
x q[82];
cx q[34], q[33];
cx q[70], q[71];
cx q[145], q[160];
cx q[83], q[68];
cx q[69], q[54];
swap q[97], q[112];
swap q[52], q[67];
swap q[128], q[143];
swap q[55], q[56];
swap q[141], q[142];
cx q[166], q[167];
cx q[196], q[181];
cx q[82], q[67];
cx q[34], q[33];
cx q[150], q[165];
x q[55];
cx q[126], q[141];
cx q[72], q[71];
cx q[143], q[158];
swap q[83], q[84];
cx q[196], q[197];
x q[82];
x q[166];
cx q[67], q[66];
x q[55];
cx q[165], q[180];
swap q[150], q[151];
swap q[152], q[167];
swap q[32], q[33];
cx q[182], q[167];
x q[152];
x q[166];
cx q[33], q[34];
cx q[81], q[82];
x q[32];
swap q[65], q[66];
swap q[67], q[68];
cx q[182], q[181];
x q[66];
cx q[153], q[152];
x q[68];
cx q[19], q[34];
cx q[31], q[32];
swap q[82], q[83];
cx q[196], q[181];
cx q[153], q[168];
cx q[152], q[137];
cx q[65], q[66];
x q[68];
cx q[19], q[20];
cx q[33], q[34];
cx q[97], q[82];
swap q[182], q[183];
x q[182];
cx q[198], q[183];
cx q[195], q[196];
cx q[153], q[168];
cx q[83], q[68];
cx q[181], q[166];
swap q[64], q[65];
swap q[136], q[137];
x q[182];
x q[183];
cx q[65], q[80];
cx q[136], q[135];
cx q[210], q[195];
x q[168];
cx q[182], q[197];
x q[80];
x q[65];
x q[136];
cx q[210], q[211];
cx q[135], q[150];
cx q[196], q[197];
cx q[65], q[64];
cx q[182], q[197];
cx q[195], q[196];
swap q[64], q[79];
x q[64];
swap q[167], q[182];
cx q[183], q[182];
cx q[79], q[64];
cx q[183], q[182];

// measurement
measure q[173]->c[0];
measure q[179]->c[1];
measure q[133]->c[2];
measure q[12]->c[3];
measure q[104]->c[4];
measure q[96]->c[5];
measure q[207]->c[6];
measure q[123]->c[7];
measure q[5]->c[8];
measure q[102]->c[9];
measure q[88]->c[10];
measure q[52]->c[11];
measure q[218]->c[12];
measure q[118]->c[13];
measure q[222]->c[14];
measure q[155]->c[15];
measure q[44]->c[16];
measure q[206]->c[17];
measure q[127]->c[18];
measure q[162]->c[19];
measure q[78]->c[20];
measure q[174]->c[21];
measure q[175]->c[22];
measure q[192]->c[23];
measure q[87]->c[24];
measure q[75]->c[25];
measure q[141]->c[26];
measure q[109]->c[27];
measure q[105]->c[28];
measure q[51]->c[29];
measure q[9]->c[30];
measure q[64]->c[31];
measure q[166]->c[32];
measure q[70]->c[33];
measure q[19]->c[34];
measure q[55]->c[35];
measure q[20]->c[36];
measure q[134]->c[37];
measure q[13]->c[38];
measure q[4]->c[39];
measure q[110]->c[40];
measure q[117]->c[41];
measure q[156]->c[42];
measure q[126]->c[43];
measure q[195]->c[44];
measure q[217]->c[45];
measure q[81]->c[46];
measure q[157]->c[47];
measure q[160]->c[48];
measure q[184]->c[49];
measure q[219]->c[50];
measure q[16]->c[51];
measure q[42]->c[52];
measure q[43]->c[53];
measure q[62]->c[54];
measure q[114]->c[55];
measure q[86]->c[56];
measure q[139]->c[57];
measure q[204]->c[58];
measure q[202]->c[59];
measure q[79]->c[60];
measure q[27]->c[61];
measure q[6]->c[62];
measure q[189]->c[63];
measure q[38]->c[64];
measure q[72]->c[65];
measure q[199]->c[66];
measure q[221]->c[67];
measure q[97]->c[68];
measure q[220]->c[69];
measure q[83]->c[70];
measure q[47]->c[71];
measure q[53]->c[72];
measure q[183]->c[73];
measure q[164]->c[74];
measure q[85]->c[75];
measure q[186]->c[76];
measure q[46]->c[77];
measure q[168]->c[78];
measure q[135]->c[79];
measure q[201]->c[80];
measure q[116]->c[81];
measure q[211]->c[82];
measure q[49]->c[83];
measure q[149]->c[84];
measure q[17]->c[85];
measure q[119]->c[86];
measure q[10]->c[87];
measure q[8]->c[88];
measure q[187]->c[89];
measure q[200]->c[90];
measure q[31]->c[91];
measure q[167]->c[92];
measure q[11]->c[93];
measure q[56]->c[94];
measure q[50]->c[95];
measure q[213]->c[96];
measure q[152]->c[97];
measure q[54]->c[98];
measure q[151]->c[99];
measure q[76]->c[100];
measure q[14]->c[101];
measure q[100]->c[102];
measure q[93]->c[103];
measure q[67]->c[104];
measure q[131]->c[105];
measure q[165]->c[106];
measure q[137]->c[107];
measure q[107]->c[108];
measure q[65]->c[109];
measure q[94]->c[110];
measure q[77]->c[111];
measure q[68]->c[112];
measure q[148]->c[113];
measure q[198]->c[114];
measure q[23]->c[115];
measure q[120]->c[116];
measure q[209]->c[117];
measure q[1]->c[118];
measure q[159]->c[119];
measure q[63]->c[120];
measure q[89]->c[121];
measure q[84]->c[122];
measure q[74]->c[123];
measure q[35]->c[124];
measure q[125]->c[125];
measure q[101]->c[126];
measure q[34]->c[127];
measure q[216]->c[128];
measure q[15]->c[129];
measure q[82]->c[130];
measure q[163]->c[131];
measure q[58]->c[132];
measure q[130]->c[133];
measure q[208]->c[134];
measure q[146]->c[135];
measure q[3]->c[136];
measure q[223]->c[137];
measure q[180]->c[138];
measure q[210]->c[139];
measure q[48]->c[140];
measure q[205]->c[141];
measure q[132]->c[142];
measure q[2]->c[143];
measure q[91]->c[144];
measure q[25]->c[145];
measure q[171]->c[146];
measure q[73]->c[147];
measure q[66]->c[148];
measure q[37]->c[149];
measure q[214]->c[150];
measure q[177]->c[151];
measure q[158]->c[152];
measure q[185]->c[153];
measure q[176]->c[154];
measure q[112]->c[155];
measure q[197]->c[156];
measure q[212]->c[157];
measure q[106]->c[158];
measure q[122]->c[159];
measure q[115]->c[160];
measure q[80]->c[161];
measure q[57]->c[162];
measure q[154]->c[163];
measure q[111]->c[164];
measure q[90]->c[165];
measure q[150]->c[166];
measure q[188]->c[167];
measure q[169]->c[168];
measure q[190]->c[169];
measure q[0]->c[170];
measure q[181]->c[171];
measure q[147]->c[172];
measure q[69]->c[173];
measure q[29]->c[174];
measure q[41]->c[175];
measure q[172]->c[176];
measure q[39]->c[177];
measure q[18]->c[178];
measure q[59]->c[179];
measure q[191]->c[180];
measure q[215]->c[181];
measure q[36]->c[182];
measure q[7]->c[183];
measure q[136]->c[184];
measure q[98]->c[185];
measure q[32]->c[186];
measure q[129]->c[187];
measure q[143]->c[188];
measure q[45]->c[189];
measure q[170]->c[190];
measure q[99]->c[191];
measure q[30]->c[192];
measure q[203]->c[193];
measure q[21]->c[194];
measure q[60]->c[195];
measure q[71]->c[196];
measure q[26]->c[197];
measure q[40]->c[198];
measure q[24]->c[199];
measure q[95]->c[200];
measure q[33]->c[201];
measure q[103]->c[202];
measure q[61]->c[203];
measure q[224]->c[204];
measure q[138]->c[205];
measure q[28]->c[206];
measure q[108]->c[207];
measure q[128]->c[208];
measure q[193]->c[209];
measure q[194]->c[210];
measure q[144]->c[211];
measure q[92]->c[212];
measure q[113]->c[213];
measure q[178]->c[214];
measure q[142]->c[215];
measure q[140]->c[216];
measure q[153]->c[217];
measure q[182]->c[218];
measure q[121]->c[219];
measure q[161]->c[220];
measure q[22]->c[221];
measure q[124]->c[222];
measure q[145]->c[223];
measure q[196]->c[224];
