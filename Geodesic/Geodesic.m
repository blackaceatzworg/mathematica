(* ::Package:: *)

(* :Mathematica Version: 8.0 *)

(* :Name: Geodesic` *)

(* :Title: Alternative Geodesic Calculation Package *)

(* :Author: Kei Misawa *)

(* :Summary:
  This package provides 3 functions, GeoDestination2, GeoDistance2 and
 GeoDirection2, which are substitute for GeoDestination, GeoDistance and
 GeoDirection respectively. This package can solve geodesic problems
 more accurately with arbitrary precision calculation by
 designating option WorkingPrecision.

  Although this package and GeographicLib refer to the same paper by
 Charles F.F.Karney[1], they differ in the way to calculate.
 GeographicLib uses functional expansion to calculate elliptic
 integrals, but this package uses Mathematica built-in functions such as 
 EllipticE. Additionally, inverse function of elliptic integral is
 represented by InverseFunction, and FindRoot is used to solve the
 inverse problems. This enables arbitrary calculation with Mathematica.

  GeoDestination2, GeoDistance2 and GeoDirection2 can be used almost
 samely as GeoDestination, GeoDistance and GeoDirection.
 However, note that not all patterns of arguments are supported.
 For example, GeoPosition and datum strings such as "WGS84" are
 not supported currently. Please use "a" and "f" options to set
 ellipsoid parameters.

  In order to calculate with arbitrary precision, set option value to
 WorkingPrecision.

 This software is released under the MIT License.
  http://opensource.org/licenses/mit-license.php
*)

(* :Context: Geodesic` *)

(* :Package Version: 0.1 *)

(* :Copyright: Copyright 2014, Kei Misawa *)

(* :Keywords: Geodesy *)

(* :Warning:  *)

(* :Sources: 
  [1] Charles F.F. Karney, "Algorithms for geodesics," J. Geodesy 87, 43-55 (2013),
      http://link.springer.com/content/pdf/10.1007%2 Fs00190-012-0578-z.pdf
*)


BeginPackage["Geodesic`"];


GeoDestination2::usage="GeoDestination2[\!\(\*
StyleBox[\"pos\", \"TI\"]\),{\!\(\*
StyleBox[\"d\", \"TI\"]\),\!\(\*
StyleBox[\"\[Alpha]\", \"TR\"]\)}] returns position from \!\(\*
StyleBox[\"pos\", \"TI\"]\) to direction \!\(\*
StyleBox[\"\[Alpha]\", \"TR\"]\) and distance \!\(\*
StyleBox[\"d\", \"TI\"]\).";


GeoDistance2::usage="GeoDistance2[{\!\(\*SubscriptBox[
StyleBox[\"lat\", \"TI\"], 
StyleBox[\"1\", \"TR\"]]\),\!\(\*SubscriptBox[
StyleBox[\"long\", \"TI\"], 
StyleBox[\"1\", \"TR\"]]\)},{\!\(\*SubscriptBox[
StyleBox[\"lat\", \"TI\"], 
StyleBox[\"2\", \"TR\"]]\),\!\(\*SubscriptBox[
StyleBox[\"long\", \"TI\"], 
StyleBox[\"2\", \"TR\"]]\)}] returns geodetic distance between 2 points.";


GeoDirection2::usage="GeoDirection2[{\!\(\*SubscriptBox[
StyleBox[\"lat\", \"TI\"], 
StyleBox[\"1\", \"TR\"]]\),\!\(\*SubscriptBox[
StyleBox[\"long\", \"TI\"], 
StyleBox[\"1\", \"TR\"]]\)},{\!\(\*SubscriptBox[
StyleBox[\"lat\", \"TI\"], 
StyleBox[\"2\", \"TR\"]]\),\!\(\*SubscriptBox[
StyleBox[\"long\", \"TI\"], 
StyleBox[\"2\", \"TR\"]]\)}] returns direction from 1st point to 2nd point.";


Begin["`Private`"];


Options[GeoDestination2]={
	"a"->6378137,
	"f"->1/298.257223563,
	WorkingPrecision->MachinePrecision
};
GeoDestination2[{lat_?NumericQ,lon_?NumericQ},{d_?NumericQ,\[Lambda]_?NumericQ},OptionsPattern[]]:=Module[{
	a=SetPrecision[OptionValue["a"],OptionValue[WorkingPrecision]],
	f=SetPrecision[OptionValue["f"],OptionValue[WorkingPrecision]],
	b,e2,eprime2,\[Phi]1,\[Alpha]1,s12,\[Beta]1,
	\[Alpha]0,\[Sigma]1,\[Omega]1,k2,\[Epsilon],s1,s2,\[Sigma]2,\[Alpha]2,
	\[Beta]2,\[Omega]2,\[Lambda]1,\[Lambda]2,\[Lambda]12,\[Phi]2
},
	b=(1-f)a;
	e2=f(2-f);
	eprime2=e2/(1-e2);
	\[Phi]1=SetPrecision[lat Degree,OptionValue[WorkingPrecision]];
	\[Alpha]1=SetPrecision[\[Lambda] Degree,OptionValue[WorkingPrecision]];
	s12=SetPrecision[d,OptionValue[WorkingPrecision]];
	\[Beta]1=ArcTan[(1-f)Tan[\[Phi]1]];
	\[Alpha]0=ArcSin[Sin[\[Alpha]1]Cos[\[Beta]1]];
	\[Sigma]1=ArcTan[Cos[\[Alpha]1]Cos[\[Beta]1],Sin[\[Beta]1]];
	\[Omega]1=ArcTan[Cos[\[Sigma]1],Sin[\[Alpha]0]Sin[\[Sigma]1]];
	k2=eprime2 Cos[\[Alpha]0]^2;
	\[Epsilon]=(Sqrt[1+k2]-1)/(Sqrt[1+k2]+1);
	s1=b EllipticE[\[Sigma]1,-k2];
	s2=s1+s12;
	\[Sigma]2=InverseFunction[EllipticE,1,2][s2/b,-k2];
	\[Alpha]2=ArcTan[Cos[\[Alpha]0]Cos[\[Sigma]2],Sin[\[Alpha]0]];
	\[Beta]2=ArcTan[Sqrt[(Cos[\[Alpha]0]Cos[\[Sigma]2])^2+Sin[\[Alpha]0]^2],Cos[\[Alpha]0]Sin[\[Sigma]2]];
	\[Omega]2=ArcTan[Cos[\[Sigma]2],Sin[\[Alpha]0]Sin[\[Sigma]2]];
	\[Lambda]12=\[Omega]2-\[Omega]1-f Sin[\[Alpha]0]I3[\[Sigma]1,\[Sigma]2,f,k2,WorkingPrecision->OptionValue[WorkingPrecision]];
	\[Phi]2=ArcTan[Tan[\[Beta]2]/(1-f)];
	{\[Phi]2/Degree,Mod[\[Lambda]12/Degree+lon,360](*,\[Alpha]2/Degree*)}
];
SyntaxInformation[GeoDestination2]={"ArgumentsPattern"->{{_,_},{_,_},OptionsPattern[]}};


(* internal function which solves inverse problems *)
Options[geoInverse]={
	"a"->6378137,
	"f"->1/298.257223563,
	WorkingPrecision->MachinePrecision,
	AccuracyGoal->Automatic,
	PrecisionGoal->Automatic,
	MaxIterations->100
};
geoInverse[{lat1_,lon1_},{lat2_,lon2_},opts:OptionsPattern[]]:=Catch[Module[{
	a=SetPrecision[OptionValue["a"],OptionValue[WorkingPrecision]],
	f=SetPrecision[OptionValue["f"],OptionValue[WorkingPrecision]],
	b,e2,eprime2,inverselat=False,inverselon=False,swapped=False,iterations=0,
	\[Phi]1,\[Phi]2,\[Lambda]12org,\[Beta]1,\[Beta]2,\[CapitalDelta],x,y,\[Mu],\[Alpha]1,\[Alpha]10,
	\[Alpha]0,\[Sigma]1,\[Omega]1,k2,\[Epsilon],\[Alpha]2,\[Sigma]2,\[Omega]2,agoal,pgoal
	s1,s2,s12
},
	b=(1-f)a;
	e2=f(2-f);
	eprime2=e2/(1-e2);
	agoal=SetPrecision[10^-If[TrueQ[OptionValue[AccuracyGoal]==Automatic],
		0.5OptionValue[WorkingPrecision]
		,
		OptionValue[AccuracyGoal]
	],OptionValue[WorkingPrecision]];
	pgoal=SetPrecision[10^-If[TrueQ[OptionValue[PrecisionGoal]==Automatic],
		0.5OptionValue[WorkingPrecision]
		,
		OptionValue[PrecisionGoal]
	],OptionValue[WorkingPrecision]];

	\[Phi]1=SetPrecision[lat1 Degree,OptionValue[WorkingPrecision]];
	\[Phi]2=SetPrecision[lat2 Degree,OptionValue[WorkingPrecision]];
	\[Lambda]12org=SetPrecision[(lon2-lon1)Degree,OptionValue[WorkingPrecision]];
	If[TrueQ[Abs[\[Phi]1]<Abs[\[Phi]2]],
		{\[Phi]1,\[Phi]2}={\[Phi]2,\[Phi]1};
		swapped=True;
	];
	If[TrueQ[\[Phi]1>0],
		\[Phi]1=-\[Phi]1;
		\[Phi]2=-\[Phi]2;
		inverselat=True;
	];
	If[TrueQ[lon2-lon1<0],
		\[Lambda]12org=-\[Lambda]12org;
		inverselon=True;
	];
	If[TrueQ[\[Lambda]12org>Pi],
		\[Lambda]12org=2Pi-\[Lambda]12org;
		inverselon=!inverselon;
	];
	\[Beta]1=ArcTan[(1-f)Tan[\[Phi]1]];
	\[Beta]2=ArcTan[(1-f)Tan[\[Phi]2]];
	(* Initial guess for \[Alpha]1 *)
	\[CapitalDelta]=(f a \[Pi] Cos[\[Beta]1]^2);
	x=(a Cos[\[Beta]1])(\[Lambda]12org-\[Pi])/\[CapitalDelta];
	y=a (\[Beta]2+\[Beta]1)/\[CapitalDelta];
	If[TrueQ[Abs[y]>agoal],
		\[Mu]=Select[Block[{\[Mu]},\[Mu]/.Solve[\[Mu]^4+2\[Mu]^3+(1-x^2-y^2)\[Mu]^2-2y^2\[Mu]-y^2==0&&\[Mu]>=0,\[Mu]]],(#>=0)&];
		If[Length[\[Mu]]>=1,
			\[Mu]=First[\[Mu]]
			,
			Throw[{$Failed,$Failed}]
		];
		\[Alpha]10=ArcTan[(y/\[Mu]),-x/(1+\[Mu])];
		,
		\[Alpha]10=ArcTan[Sqrt[Max[0,1-x^2]],-x];
	];

	(* Solve \[Alpha]1 *)
	\[Alpha]1=\[Alpha]1/.First@FindRoot[
		-ArcTan[Cos[ArcTan[Cos[\[Alpha]1] Cos[\[Beta]1],Sin[\[Beta]1]]],Cos[\[Beta]1] Sin[\[Alpha]1] Sin[ArcTan[Cos[\[Alpha]1] Cos[\[Beta]1],Sin[\[Beta]1]]]]
		+ArcTan[Cos[ArcTan[Sqrt[-Cos[\[Beta]1]^2+Cos[\[Alpha]1]^2 Cos[\[Beta]1]^2+Cos[\[Beta]2]^2],Sin[\[Beta]2]]],Cos[\[Beta]1] Sin[\[Alpha]1] Sin[ArcTan[Sqrt[-Cos[\[Beta]1]^2+Cos[\[Alpha]1]^2 Cos[\[Beta]1]^2+Cos[\[Beta]2]^2],Sin[\[Beta]2]]]]
		-f Cos[\[Beta]1] I3[ArcTan[Cos[\[Alpha]1] Cos[\[Beta]1],Sin[\[Beta]1]],ArcTan[Sqrt[-Cos[\[Beta]1]^2+Cos[\[Alpha]1]^2 Cos[\[Beta]1]^2+Cos[\[Beta]2]^2],Sin[\[Beta]2]],f,eprime2 (1-Cos[\[Beta]1]^2 Sin[\[Alpha]1]^2),WorkingPrecision->OptionValue[WorkingPrecision]] Sin[\[Alpha]1]
		==\[Lambda]12org
		,{\[Alpha]1,\[Alpha]10,0,Pi},
		Evaluate@FilterRules[{
			Evaluated->False,
			Method->"Secant",
			opts},Options[FindRoot]]
	];

	\[Alpha]0=ArcSin[Sin[\[Alpha]1]Cos[\[Beta]1]];
	\[Sigma]1=ArcTan[Cos[\[Alpha]1]Cos[\[Beta]1],Sin[\[Beta]1]];
	\[Sigma]2=ArcTan[Sqrt[Cos[\[Alpha]1]^2Cos[\[Beta]1]^2+(Cos[\[Beta]2]^2-Cos[\[Beta]1]^2)],Sin[\[Beta]2]];
	\[Alpha]2=ArcTan[Cos[\[Alpha]0]Cos[\[Sigma]2],Sin[\[Alpha]0]];
	k2=eprime2 Cos[\[Alpha]0]^2;
	\[Epsilon]=(Sqrt[1+k2]-1)/(Sqrt[1+k2]+1);
	s1=Re[b EllipticE[\[Sigma]1,-k2]];
	s2=Re[b EllipticE[\[Sigma]2,-k2]];
	s12=s2-s1;

	If[swapped,{\[Alpha]1,\[Alpha]2}={Pi-\[Alpha]2,Pi-\[Alpha]1};];
	If[inverselat,{\[Alpha]1,\[Alpha]2}={Pi-\[Alpha]1,Pi-\[Alpha]2};];
	If[inverselon,{\[Alpha]1,\[Alpha]2}={-\[Alpha]1,-\[Alpha]2}];
	\[Alpha]1=Mod[\[Alpha]1,2Pi];
	\[Alpha]2=Mod[\[Alpha]2,2Pi];
	Throw[{s12,\[Alpha]1/Degree,\[Alpha]2/Degree}];
]];


(* internal function which calculates elliptic integral numerically *)
Clear[I3];
Options[I3]=Options[NIntegrate];
I3[\[Sigma]1_?NumericQ,\[Sigma]2_?NumericQ,f_,k2_,opts:OptionsPattern[]]:=NIntegrate[
	SetPrecision[(2-f)/(1+(1-f)Sqrt[1+k2 Sin[\[Sigma]0]^2]),OptionValue[WorkingPrecision]]
	,{\[Sigma]0,\[Sigma]1,\[Sigma]2},opts
]


Options[GeoDistance2]=Options[geoInverse];
GeoDistance2[{lat1_?NumericQ,lon1_?NumericQ},{lat2_?NumericQ,lon2_?NumericQ},opts:OptionsPattern[]]:=geoInverse[
{lat1,lon1},{lat2,lon2},opts][[1]];
SyntaxInformation[GeoDistance2]={"ArgumentsPattern"->{{_,_},{_,_},OptionsPattern[]}};


Options[GeoDirection2]=Options[geoInverse];
GeoDirection2[{lat1_?NumericQ,lon1_?NumericQ},{lat2_?NumericQ,lon2_?NumericQ},opts:OptionsPattern[]]:=geoInverse[
{lat1,lon1},{lat2,lon2},opts][[2]];
SyntaxInformation[GeoDirection2]={"ArgumentsPattern"->{{_,_},{_,_},OptionsPattern[]}};


End[];


EndPackage[];
