(* ::Package:: *)

(* ::Title:: *)
(*Kerr-Newman metric*)


AbsoluteTiming[

ds2Matrix[dx_,ds2_]:=Table[D[ds2,dx[[i1]],dx[[i2]]]/2,{i1,Length[dx]},{i2,Length[dx]}];
dcoo:=Table[ToExpression[StringJoin["d",ToString[coo[[i1]]]]],{i1,Length[coo]}];
(**)
gxx[i1_,i2_]:=metric[[i1,i2]];
gXX[i1_,i2_]:=Inverse[metric][[i1,i2]];
GXxx[i1_,i2_,i3_]:=Sum[gXX[i1,ia](-D[gxx[i2,i3],coo[[ia]]]+D[gxx[i3,ia],coo[[i2]]]+D[gxx[ia,i2],coo[[i3]]]),{ia,dim}]/2;
RxxxX[i1_,i2_,i3_,i4_]:=D[GXxx[i4,i1,i3],coo[[i2]]]-D[GXxx[i4,i2,i3],coo[[i1]]]+Sum[GXxx[ia,i3,i1]GXxx[i4,i2,ia]-GXxx[ia,i3,i2]GXxx[i4,i1,ia],{ia,dim}];
Rxx[i1_,i3_]:=Sum[RxxxX[i1,ia,i3,ia],{ia,dim}];
R[]:=Sum[gXX[ia,ib]Rxx[ia,ib],{ia,dim},{ib,dim}];
(**)
U1T1[T_][i1_]:=Sum[gXX[i1,ia]T[ia],{ia,dim}];
U1T2[T_][i1_,i2_]:=Sum[gXX[i1,ia]T[ia,i2],{ia,dim}];
U2T2[T_][i1_,i2_]:=Sum[gXX[i2,ia]T[i1,ia],{ia,dim}];
U1T3[T_][i1_,i2_,i3_]:=Sum[gXX[i1,ia]T[ia,i2,i3],{ia,dim}];
U2T3[T_][i1_,i2_,i3_]:=Sum[gXX[i2,ia]T[i1,ia,i3],{ia,dim}];
U3T3[T_][i1_,i2_,i3_]:=Sum[gXX[i3,ia]T[i1,i2,ia],{ia,dim}];
U1T4[T_][i1_,i2_,i3_,i4_]:=Sum[gXX[i1,ia]T[ia,i2,i3,i4],{ia,dim}];
U2T4[T_][i1_,i2_,i3_,i4_]:=Sum[gXX[i2,ia]T[i1,ia,i3,i4],{ia,dim}];
U3T4[T_][i1_,i2_,i3_,i4_]:=Sum[gXX[i3,ia]T[i1,i2,ia,i4],{ia,dim}];
U4T4[T_][i1_,i2_,i3_,i4_]:=Sum[gXX[i4,ia]T[i1,i2,i3,ia],{ia,dim}];
(**)
D1T1[T_][i1_]:=Sum[gxx[i1,ia]T[ia],{ia,dim}];
D1T2[T_][i1_,i2_]:=Sum[gxx[i1,ia]T[ia,i2],{ia,dim}];
D2T2[T_][i1_,i2_]:=Sum[gxx[i2,ia]T[i1,ia],{ia,dim}];
D1T3[T_][i1_,i2_,i3_]:=Sum[gxx[i1,ia]T[ia,i2,i3],{ia,dim}];
D2T3[T_][i1_,i2_,i3_]:=Sum[gxx[i2,ia]T[i1,ia,i3],{ia,dim}];
D3T3[T_][i1_,i2_,i3_]:=Sum[gxx[i3,ia]T[i1,i2,ia],{ia,dim}];
D1T4[T_][i1_,i2_,i3_,i4_]:=Sum[gxx[i1,ia]T[ia,i2,i3,i4],{ia,dim}];
D2T4[T_][i1_,i2_,i3_,i4_]:=Sum[gxx[i2,ia]T[i1,ia,i3,i4],{ia,dim}];
D3T4[T_][i1_,i2_,i3_,i4_]:=Sum[gxx[i3,ia]T[i1,i2,ia,i4],{ia,dim}];
D4T4[T_][i1_,i2_,i3_,i4_]:=Sum[gxx[i4,ia]T[i1,i2,i3,ia],{ia,dim}];
(**)
DxT[T_][i1_]:=D[T,coo[[i1]]];
DxTx[Tx_][i1_,i2_]:=D[Tx[i2],coo[[i1]]]+Sum[-GXxx[ia,i1,i2]Tx[ia],{ia,dim}];
DxTX[TX_][i1_,i2_]:=D[TX[i2],coo[[i1]]]+Sum[+GXxx[i2,i1,ia]TX[ia],{ia,dim}];
DxTxx[Txx_][i1_,i2_,i3_]:=D[Txx[i2,i3],coo[[i1]]]+Sum[-GXxx[ia,i1,i2]Txx[ia,i3]-GXxx[ia,i1,i3]Txx[i2,ia],{ia,dim}];
DxTXx[TXx_][i1_,i2_,i3_]:=D[TXx[i2,i3],coo[[i1]]]+Sum[+GXxx[i2,i1,ia]TXx[ia,i3]-GXxx[ia,i1,i3]TXx[i2,ia],{ia,dim}];
DxTxX[TxX_][i1_,i2_,i3_]:=D[TxX[i2,i3],coo[[i1]]]+Sum[-GXxx[ia,i1,i2]TxX[ia,i3]+GXxx[i3,i1,ia]TxX[i2,ia],{ia,dim}];
DxTXX[TXX_][i1_,i2_,i3_]:=D[TXX[i2,i3],coo[[i1]]]+Sum[+GXxx[i2,i1,ia]TXX[ia,i3]+GXxx[i3,i1,ia]TXX[i2,ia],{ia,dim}];
DxTxxx[Txxx_][i1_,i2_,i3_,i4_]:=D[Txxx[i2,i3,i4],coo[[i1]]]+Sum[-GXxx[ia,i1,i2]Txxx[ia,i3,i4]-GXxx[ia,i1,i3]Txxx[i2,ia,i4]-GXxx[ia,i1,i4]Txxx[i2,i3,ia],{ia,dim}];
DxTXxx[TXxx_][i1_,i2_,i3_,i4_]:=D[TXxx[i2,i3,i4],coo[[i1]]]+Sum[+GXxx[i2,i1,ia]TXxx[ia,i3,i4]-GXxx[ia,i1,i3]TXxx[i2,ia,i4]-GXxx[ia,i1,i4]TXxx[i2,i3,ia],{ia,dim}];
DxTxXx[TxXx_][i1_,i2_,i3_,i4_]:=D[TxXx[i2,i3,i4],coo[[i1]]]+Sum[-GXxx[ia,i1,i2]TxXx[ia,i3,i4]+GXxx[i3,i1,ia]TxXx[i2,ia,i4]-GXxx[ia,i1,i4]TxXx[i2,i3,ia],{ia,dim}];
DxTxxX[TxxX_][i1_,i2_,i3_,i4_]:=D[TxxX[i2,i3,i4],coo[[i1]]]+Sum[-GXxx[ia,i1,i2]TxxX[ia,i3,i4]-GXxx[ia,i1,i3]TxxX[i2,ia,i4]+GXxx[i4,i1,ia]TxxX[i2,i3,ia],{ia,dim}];
DxTXXx[TXXx_][i1_,i2_,i3_,i4_]:=D[TXXx[i2,i3,i4],coo[[i1]]]+Sum[+GXxx[i2,i1,ia]TXXx[ia,i3,i4]+GXxx[i3,i1,ia]TXXx[i2,ia,i4]-GXxx[ia,i1,i4]TXXx[i2,i3,ia],{ia,dim}];
DxTXxX[TXxX_][i1_,i2_,i3_,i4_]:=D[TXxX[i2,i3,i4],coo[[i1]]]+Sum[+GXxx[i2,i1,ia]TXxX[ia,i3,i4]-GXxx[ia,i1,i3]TXxX[i2,ia,i4]+GXxx[i4,i1,ia]TXxX[i2,i3,ia],{ia,dim}];
DxTxXX[TxXX_][i1_,i2_,i3_,i4_]:=D[TxXX[i2,i3,i4],coo[[i1]]]+Sum[-GXxx[ia,i1,i2]TxXX[ia,i3,i4]+GXxx[i3,i1,ia]TxXX[i2,ia,i4]+GXxx[i4,i1,ia]TxXX[i2,i3,ia],{ia,dim}];
DxTXXX[TXXX_][i1_,i2_,i3_,i4_]:=D[TXXX[i2,i3,i4],coo[[i1]]]+Sum[+GXxx[i2,i1,ia]TXXX[ia,i3,i4]+GXxx[i3,i1,ia]TXXX[i2,ia,i4]+GXxx[i4,i1,ia]TXXX[i2,i3,ia],{ia,dim}];
(**)
DXT[T_]:=U1T1[DxT[T]];
DXTx[T_]:=U1T2[DxTx[T]];
DXTX[T_]:=U1T2[DxTX[T]];
DXTxx[T_]:=U1T3[DxTxx[T]];
DXTXx[T_]:=U1T3[DxTXx[T]];
DXTxX[T_]:=U1T3[DxTxX[T]];
DXTXX[T_]:=U1T3[DxTXX[T]];
DXTxxx[T_]:=U1T4[DxTxxx[T]];
DXTXxx[T_]:=U1T4[DxTXxx[T]];
DXTxXx[T_]:=U1T4[DxTxXx[T]];
DXTxxX[T_]:=U1T4[DxTxxX[T]];
DXTXXx[T_]:=U1T4[DxTXXx[T]];
DXTXxX[T_]:=U1T4[DxTXxX[T]];
DXTxXX[T_]:=U1T4[DxTxXX[T]];
DXTXXX[T_]:=U1T4[DxTXXX[T]];
DxDxT[T_][i1_,i2_]:=DxTx[DxT[T]][i1,i2];
boxT[T_]:=Sum[U1T2[DxDxT[T]][ia,ia],{ia,dim}];
(**)
print1[T_]:=Do[{{i1},T[i1]}//Simplify//Print,{i1,dim}];
print2[T_]:=Do[{{i1,i2},T[i1,i2]}//Simplify//Print,{i1,dim},{i2,dim}];
print3[T_]:=Do[{{i1,i2,i3},T[i1,i2,i3]}//Simplify//Print,{i1,dim},{i2,dim},{i3,dim}];
print4[T_]:=Do[{{i1,i2,i3,i4},T[i1,i2,i3,i4]}//Simplify//Print,{i1,dim},{i2,dim},{i3,dim},{i4,dim}];
matrix[T2_]:=Table[T2[i1,i2],{i1,dim},{i2,dim}]//MatrixForm;
(******************************************************************************************************************************************************************)
coo={t,r,\[Theta],\[Phi]};
dcoo=dcoo;
dim=Length[coo];
ds2=-(1-(2M r-Q^2)/\[Rho][r,\[Theta]]^2)dt^2+\[Rho][r,\[Theta]]^2/\[CapitalDelta][r]dr^2+\[Rho][r,\[Theta]]^2d\[Theta]^2+((r^2+a^2)^2-\[CapitalDelta][r]a^2Sin[\[Theta]]^2)Sin[\[Theta]]^2/\[Rho][r,\[Theta]]^2d\[Phi]^2-2a(2M r-Q^2)Sin[\[Theta]]^2/\[Rho][r,\[Theta]]^2dt d\[Phi]; 
\[Rho][r,\[Theta]]=Sqrt[r^2+a^2Cos[\[Theta]]^2]; 
\[CapitalDelta][r]=r^2-2M r+a^2+Q^2; 
metric=ds2Matrix[dcoo,ds2]//Simplify;
Do[gxx[i1,i2]=gxx[i1,i2]//Simplify,{i1,dim},{i2,dim}];
Do[gXX[i1,i2]=gXX[i1,i2]//Simplify,{i1,dim},{i2,dim}];
tGXxx=Table[Print[{i1,i2,i3}];GXxx[i1,i2,i3]//Simplify[#,TimeConstraint->Infinity]&,{i1,dim},{i2,dim},{i3,dim}]//Parallelize;
Do[GXxx[i1,i2,i3]=tGXxx[[i1,i2,i3]],{i1,dim},{i2,dim},{i3,dim}];
tRxxxX=Table[Print[{i1,i2,i3,i4}];RxxxX[i1,i2,i3,i4]//Simplify[#,TimeConstraint->Infinity]&,{i1,dim},{i2,dim},{i3,dim},{i4,dim}]//Parallelize;
Do[RxxxX[i1,i2,i3,i4]=tRxxxX[[i1,i2,i3,i4]],{i1,dim},{i2,dim},{i3,dim},{i4,dim}];
tRxx=Table[Print[{i1,i2}];Rxx[i1,i2]//Simplify[#,TimeConstraint->Infinity]&,{i1,dim},{i2,dim}]//Parallelize;
Do[Rxx[i1,i2]=tRxx[[i1,i2]],{i1,dim},{i2,dim}];
R[]=R[]//Simplify[#,TimeConstraint->Infinity]&

]
(*i5-6600,        MMA12.1.1.0, 20s, 2252s*)
(*MacbookPro2015, MMA11.0.1.0, 25s, 3074s*)