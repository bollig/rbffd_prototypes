(* ::Package:: *)

(**************************************************************************
*                                                                         *
*                               Ricci                                     *
*        A Mathematica package for doing tensor calculations              *
*                      in differential geometry                           *
*                                                                         *
*                            Version 1.53                                 *
*                                                                         *
*                                                                         *
*                           By John M. Lee                                *
*    assisted by Dale Lear, John Roth, Lee Nave, and Larry Peterson       *
*                                                                         *
*                Copyright (c) 1992 - 2011 John M. Lee                    *
*                        All rights reserved                              *
*                                                                         *
*         Development of this software was supported in part              *
*                by NSF grants DMS-9101832, DMS-9404107                   *
*                                                                         *
*                                                                         *
* This software package and its accompanying documentation are provided   *
* as is, without guarantee of support or maintenance.  The copyright      *
* holder makes no express or implied warranty of any kind with respect    *
* to this software, including implied warranties of merchantability or    *
* fitness for a particular purpose, and is not liable for any damages     *
* resulting in any way from its use.                                      *
*                                                                         *
* Everyone is granted permission to copy, modify and redistribute this    *
* software package and its accompanying documentation, provided that:     *
*  1. All copies contain this notice in the main program file and in the  *
*     supporting documentation.                                           *
*  2. All modified copies carry a prominent notice stating who made the   *
*     last modification and the date of such modification.                *
*  3. No charge is made for this software or works derived from it, with  *
*     the exception of a distribution fee to cover the cost of materials  *
*     and/or transmission.                                                *
*                                                                         *
**************************************************************************)

BeginPackage["Ricci`"];

Print[" -- Ricci Version 1.53 (September 16, 2011) --"];
Print["    Copyright 1992 - 2011 John M. Lee"];
Print["    Problem reports or suggestions to:"];
Print["      lee@math.washington.edu"];

If [ $VersionNumber >= 3.0,
     Print[" -- For properly formatted output in a notebook,"];
     Print["    set the default output format type to OutputForm"]
    ];


(* This function is used to append our messages to the system messages
   for system functions that are modified by Ricci.  If the system
   message doesn't have head String, then ignore it.  Otherwise, see
   if our message is already there.  If so, leave it alone; if not,
   append our message to the existing one. *)

Begin["`Private`"];

AppendMessage[ sysmsg_, string_ ] :=
  Which[
    Head[sysmsg] =!= String,
      string,
    StringMatchQ[ sysmsg, Evaluate[ "*" <> string <> "*" ] ],
      sysmsg,
    True,
      sysmsg <> "\n\n" <> string];

End[ (* "Ricci`Private`" *) ];

(*****************************************************************)

(** Initialize usage messages.  These messages define all the
    user-accessible symbols, defined in the "Ricci`" context. **)

AbsorbMetrics::usage = "AbsorbMetrics[expr] simplifies expr by
eliminating any metrics that are contracted with other tensors, and
using them to raise or lower indices. AbsorbMetrics[x,n] or
AbsorbMetrics[x,{n1,n2,...}] applies AbsorbMetrics only to term n or
to terms n1,n2,... of x.
\nOption:
\n* Mode -> All or NoOneDims. NoOneDims means that AbsorbMetrics
should not absorb one-dimensional metrics (unless they are paired with
other metrics). Default is All.";

Alt::usage = "Alt[x] is the alternating part of the tensor expression
x.";

Alternating::usage = "Alternating is a value for the Symmetries option
of DefineTensor.";

Any::usage = "Any can be used as a value for the Bundle option of the
DefineTensor command. Bundle -> Any means that indices from any
bundle can be inserted.";

Bar::usage =
"The internal form for a barred index is L[i[Bar]] or U[i[Bar]]. In
input form, these can be abbreviated LB[i] and UB[i].
\nThe
internal form for the conjugate of a tensor is Tensor[name[Bar],{...},{...}].
In input form, this is typed Conjugate[name] [...] [...].";

Basis::usage = "Basis is the name used for generic basis vectors and
covectors for any bundle. Basis vectors are generated automatically
by BasisExpand. For example, if index i is associated with bundle b,
then Basis[L[i]] and Basis[U[i]] represent contravariant and covariant
basis elements for b (i.e., basis elements for b and its dual),
respectively; Basis[LB[i]] and Basis[UB[i]] represent basis elements
for the conjugate bundle Conjugate[b].
\nIf you insert indices into
Basis vectors, you get Kronecker delta functions, metrics, or zero as
appropriate.";

BasisExpand::usage = "BasisExpand[x] converts x to a sum of component
expressions multiplied by basis vectors and covectors, Basis[L[i]] and
Basis[U[j]]. BasisExpand[x,n] or BasisExpand[x,{n1,n2,...}] applies
BasisExpand only to term n or terms n1,n2,... of x.";

BasisGather::usage = "BasisGather[x,tensor] attempts to recognize the
basis expression for \"tensor\" in x, and replaces it by \"tensor\".
BasisGather[x,{tensor1,tensor2,...}] does the same thing for several
tensors at once.";

BasisName::usage = "BasisName -> name is an option for DefineBundle.
It specifies the name to be used in the input and output forms of the
basis vectors for the bundle being defined. Default = Basis.";

BianchiRules::usage = "BianchiRules[i,j,k] converts tensors with
Symmetries -> RiemannSymmetries containing the indices i,j,k in order
to two-term sums, using the first and second Bianchi identities.";

Bundle::usage = "Bundle -> name is a DefineTensor option, specifying
the name of the bundle the tensor is to be associated with. This can
be a list of bundle names, meaning that the tensor is associated with
the direct sum of the bundles in the list.
\nThe function Bundle[i]
returns the name of the bundle associated with the index i.";

BundleDummy::usage =
"BundleDummy[bundle] returns the symbol that is used for
computer-generated dummy indices associated with the bundle. This is the last
index name in the list BundleIndices[bundle]. For most bundles,
new dummy indices will be of the form \"kn\", where k=BundleDummy[bundle] and
n is an integer. For one-dimensional bundles, n is omitted.";

BundleIndices::usage =
"BundleIndices[bundle] returns a list of the index names currently associated
with bundle.";

BundleQ::usage = "BundleQ[x] returns True if x is the name of a bundle,
and False otherwise.";

Bundles::usage = "Bundles[x] returns a list of the bundles
associated with index positions in the tensor expression x; each
entry in the list is a sublist giving the allowed bundles for the
corresponding index position.";

Co::usage = "Co is an abbreviation for Covariant.";

CoBasisName::usage = "CoBasisName -> name is an option for
DefineBundle. It specifies the name to be used in input and output
form for the basis covectors for the bundle being defined. Default is
the name given in the BasisName option, or Basis if BasisName is not
specified.";

CollectConstants::usage = "CollectConstants[x] groups together terms
in the tensor expression x having the same tensor factors but
different constant factors. CollectConstants[x,n] or
CollectConstants[x,{n1,n2,...}] applies CollectConstants only to term
n or to terms n1,n2,... of x.";

CommuteCovD::usage = "CommuteCovD[ x, L[i], L[j] ] changes all
adjacent occurrences of indices L[i], L[j] after the \";\" to L[j],
L[i] by adding appropriate curvature and torsion terms.";

CommutingFrame::usage = "CommutingFrame is an option for DefineBundle,
which is meaningful only for tangent bundles and their subbundles. If
True, it means that the default frame for that bundle is always
assumed to consist of commuting vector fields.";

CompatibilityRule::usage = "CompatibilityRule is a rule that
transforms derivatives of metric components into connection
coefficients, using the fact that the default connection is compatible
with the metric.";

Complex::usage = Ricci`Private`AppendMessage[Complex::usage, "Ricci
recognizes Complex as a value for the Type option of DefineBundle,
DefineConstant, DefineMathFunction and DefineTensor."];

Con::usage = "Con is an abbreviation for Contravariant.";

Conjugate::usage = Ricci`Private`AppendMessage[Conjugate::usage, "The
Ricci package modifies Conjugate to handle tensor expressions and
indices. The behavior of a tensor, constant, index, mathematical
function, or bundle under conjugation is determined by the Type option
when the object is defined."];

Conn::usage = "Conn is the generic connection form for the default
connection on any bundle. It is defined as a product tensor of rank
{2,1}; when the first two indices are inserted, it represents the
matrix of connection 1-forms associated with the default basis. When
all three indices are inserted, it represents the Christoffel symbols
of the default connection relative to the default basis.";

Connection::usage = "Connection -> cn is an option for some of the
Ricci differentiation functions. It specifies that covariant
derivatives are to be taken with respect to the connection cn instead
of the default connection. Ordinarily, cn will be an expression of
the form \"Conn + diff\", where diff is a 3-tensor expression
representing the difference tensor between the default connection and
cn.";

ConnToMetricRule::usage = "ConnToMetricRule is a rule that causes
components of the default connection Conn to be converted to
expressions involving directional derivatives of the metric. This
rule applies only to fully-indexed connection components whose indices
all refer to a single Riemannian tangent bundle with the option
CommutingFrame -> True.";

ConstantFactor::usage = "ConstantFactor[x] returns the product
of all the constant factors in x, which should be a monomial.";

ConstantQ::usage = "ConstantQ[x] returns True if there are no explicit
tensors in x, and False otherwise.";

ContractedBianchiRules::usage = "ContractedBianchiRule is a rule that
simplifies contracted covariant derivatives of the Riemannian and
Ricci curvature tensors using contracted versions of the second
Bianchi identity.";

Contravariant::usage = "Contravariant is a value for the Variance option
of DefineTensor, and may be abbreviated Con. If an index slot is
Contravariant, indices in that slot are upper by default.";

CorrectAllVariances::usage = "CorrectAllVariances[x] changes the
variance (upper to lower or lower to upper) of indices in x whose
variances are not correct for their positions, by inserting
appropriate metric coefficients.
\nOption:
\n* Mode -> All or OneDims. OneDims means that CorrectAllVariances
should correct variances only of one-dimensional indices and indices
that appear inside differential operators such as Del or Extd.
Default is All.";

CovD::usage = "If x is a component expression (no unfilled index
slots), then x[ L[i], L[j] ] or CovD[ x, {L[i],L[j]} ] is the
component of the covariant derivative of x in the L[i],L[j]
directions. In output form, covariant derivatives of a tensor are
represented by indices following a semicolon."<>
"\nOptions (for the second format only):
\n* Connection -> cn specifies that covariant derivatives are to be
taken with respect to the connection cn instead of the default
connection.
\n* Metric -> g: a symmetric 2-tensor expression, indicating that
the covariant derivatives are to be taken with respect to the
Levi-Civita connection of g instead of the default metric for x's bundle
(which is assumed to be Riemannian).";

CovDExpand::usage = "CovDExpand[x] converts all covariant derivatives
in x to ordinary directional derivatives (represented as
Del[Basis[L[i]],...]) and connection coefficients. CovDExpand[x,n] or
CovDExpand[x,{n1,n2,...}] applies CovDExpand only to term n or to terms
n1,n2,... of x.";

CovDSimplify::usage = "CovDSimplify[x] attempts to simplify x as much
as possible by ordering all dummy indices, including those that occur
after the \";\" in component expressions. CovDSimplify[x,n] or
CovDSimplify[x,{n1,n2,...}] applies CovDSimplify only to term n or to
terms n1,n2,... of x. CovDSimplify first calls TensorSimplify, then
OrderCovD, then TensorSimplify again.";

Covariant::usage = "Covariant is a tensor Variance. If an index slot
is Covariant, indices in that slot are lower by default. May be
abbreviated Co.";

Curv::usage = "Curv is the generic curvature tensor associated with
the default connection on any bundle. It is generated when covariant
derivatives are commuted, and when SecondStructureRule is applied to
derivatives of connection forms. Curv is a product tensor with rank
{2,2}, which is Alternating in its last two indices. Inserting the
first two indices yields the matrix of curvature 2-forms associated
with the default basis. Inserting all four indices yields the
coefficients of the curvature tensor.";

Curvature::usage = "Curvature[cn] is the curvature tensor associated
to the connection cn. Ordinarily, cn will be an expression of the
form \"Conn + diff\", where diff is a 3-tensor expression representing
the difference tensor between the default connection and cn.";

CurvToConnRule::usage = "CurvToConnRule is a rule that converts
components of curvature tensors to connection coefficients and their
directional derivatives.";

Declare::usage = "Declare[name,options] can be used to change certain
options for a previously-defined bundle, constant, index, tensor, or
mathematical function. Declare[{name1,name2,...},options] changes
options for several names at once. For constants and math functions,
only the Type option is allowed. For indices, the only allowable
option is TeXFormat. For bundles, the allowable options are
FlatConnection, ParallelFrame, OrthonormalFrame, PositiveDefinite,
CommutingFrame, and TorsionFree. For tensors, the allowable options
are Type, TeXFormat, and Bundle.";

DeclareBundle::usage = Declare::usage;

DeclareConstant::usage = Declare::usage;

DeclareIndex::usage = Declare::usage;

DeclareTensor::usage = Declare::usage;

DeclareMathFunction::usage = Declare::usage;

DefineBundle::usage =
"DefineBundle[ name, dim, metric, {indices} ] defines a bundle.
\n* name: A symbol that uniquely identifies the bundle.
\n* dim: The dimension of the bundle. A positive integer or
symbolic constant.
\n* metric: A name for the bundle's metric.
\n* indices: A list of index names to be associated with the
bundle."<>
"\nOptions:
\n* Type -> Complex or Real. Default is Real.
\n* TangentBundle -> bundle. The name of the underlying tangent
bundle for the bundle being defined. Default is $DefaultTangentBundle
if defined, otherwise this bundle itself (and its conjugate if the
bundle is complex)."<> "\n* FlatConnection -> True or False. Default
is False.
\n* ParallelFrame -> True or False. Default is False.
\n* OrthonormalFrame -> True or False. Default is False.
\n* CommutingFrame -> True or False. Default is False.
\n* PositiveDefinite -> True or False. Default is True.
\n* TorsionFree -> True or False. Default is True.
\n* BasisName -> name. A name for the default basis for this bundle.
\n* CoBasisName -> name. A name for the default covariant basis for
this bundle.
\n* MetricType -> Riemannian or Normal. Default is Normal. If Riemannian
is specified, there are additional options for curvature
conventions.
\n* Quiet -> True or False. Default is False, or the value of $Quiet if set.";

DefineConstant::usage = "DefineConstant[symbol] defines symbol to be a
constant with respect to differentiation in all space variables.
\nOptions: 
\n* Type -> type. This can be a single keyword or a list of
keywords, chosen from among Real, Complex, Imaginary, Positive,
Negative, NonPositive, NonNegative, Integer, Even, or Odd. Default is
Real.
\n* Quiet -> True or False. Default is False, or the value of
$Quiet if set.";

DefineIndex::usage = "DefineIndex[ {i,j,k}, bundle ]
causes the index names i,j,k to be associated with bundle.
\n* The first argument must be a symbol or list of symbols.
\n* The second argument must be a the name of a bundle.
\nOptions:
\n* TeXFormat -> \"texformat\". Default
is the index name itself.
\n* Quiet -> True or False. Default is False, or the value of $Quiet if set.";

DefineMathFunction::usage = "DefineMathFunction[f] declares f to be a
scalar-valued function of one real or complex variable that will be
used in tensor expressions.
\nOptions:
\n* Type -> type. This can be a single keyword or a list of keywords,
chosen from among Real, Complex, Imaginary, Automatic, Positive,
Negative, NonPositive, or NonNegative. Default is Real.
\n* Quiet -> True or False. Default is False, or the value of
$Quiet if set.";

DefineRelation::usage = "DefineRelation[ tensor, expr,
cond ] defines a relation of the form
\n	tensor := expr /; cond
\n associated with the tensor\'s name. The \"cond\" is optional. The
first argument \"tensor\" represents a tensor name with or without
indices, and can be written in Ricci's input form. If index names used
in \"tensor\" are associated with bundles, then the relation will only
be applied when indices from those bundles appear in those positions.
The relation will substitute x (suitably modified) in place of
any expression that matches tensor after insertion of indices,
covariant differentiation, conjugation, or raising or lowering of certain
indices."<>
"\nOption:
\nQuiet -> True or False. Default is False, or the value of $Quiet if set.";

DefineRule::usage = "DefineRule[ rulename, lhs, rhs, cond ] defines a
rule named \"rulename\" of the form lhs :> rhs /; cond. The \"cond\"
is optional. The rule is appended to previous rules defined with the
same name, unless the option NewRule -> True is specified. If index
names used on the left-hand side are associated with bundles, then the
rule is applied only when indices from those bundles appear in
those positions. If lhs is a single tensor with or without indices,
then the rule will substitute rhs (suitably modified) in place of any
expression that matches lhs after insertion of indices, covariant
differentiation, conjugation, or raising or lowering of certain
indices."<>
"\nOptions:
\nQuiet -> True or False. Default is False or the value of $Quiet if
set.
\nNewRule -> True or False. Default is False.";

DefineTensor::usage =
"DefineTensor[ name, rank ] defines a tensor.
\n* name: A symbol that uniquely identifies the tensor.
\n* rank: Rank of the tensor. For ordinary tensors, this
must be a non-negative integer. For product tensors, this
is a list of non-negative integers." <>
"\nOptions:
\n* Symmetries -> sym:
Symmetries of the tensor.
Default is NoSymmetries.
\n* Type -> type. For rank-0 tensors, this can be a single keyword or
a list of keywords, chosen from among Complex, Real, Imaginary,
Positive, Negative, NonPositive, or NonNegative. For higher-rank
tensors, it can be Real, Complex, or Imaginary. The default is always
Real.
\n* TeXFormat -> \"texformat\":    Default is the tensor's name." <>
"\n* Bundle -> bundle:  The bundle with
which the tensor is associated. Can be a list, meaning
the tensor is associated with the direct sum of all the
bundles in the list. Can be Any, to indicate that the tensor
will accept indices from any bundle.
Default is $DefaultTangentBundle if
set; otherwise Bundle must be specified.
\n* Variance -> Covariant or Contravariant:  Default is Covariant.
This can be a list whose length is equal to the total rank
of the tensor, in which case each entry in the list specifies
the variance of the corresponding index slot. May be
abbreviated Co and Con.
\n* Quiet -> True or False:  Default is False, or the value
of $Quiet if set."

DefineTensorSymmetries::usage =
"DefineTensorSymmetries[name,{perm1,sgn1,...,permN,sgnN}] defines a
TensorSymmetry object that can be used in DefineTensor. The permi
must be lists that are nontrivial permutations of {1,2,...,d} and the
sgni must be constants (ordinarily plus or minus 1).";

Del::usage = "Del[x] is the total covariant derivative of
the tensor expression x. If x is a k-tensor, then Del[x] is a
k+1-tensor. The last index of Del[x] is the one generated by
differentiation.
\nDel[v,x] is the covariant derivative of x in
the direction v; v must be a 1-tensor expression.
\nOptions:
\n* Connection -> cn specifies that the covariant derivative is to be
taken with respect to the connection cn instead of the default
connection.
\n* Metric -> g: a symmetric 2-tensor expression, indicating that
the covariant derivative is to be taken with respect to the
Levi-Civita connection of g instead of the default metric for x's bundle
(which is assumed to be Riemannian).";

Det::usage = Ricci`Private`AppendMessage[ Det::usage,
"If x is a 2-tensor expression, Ricci interprets Det[x] as the determinant
of x, using the metric to convert x to a {Contravariant,Covariant} tensor
if necessary." ];

Dimension::usage =
"Dimension[bundle] returns the bundle\'s dimension.";

Div::usage = "Div[x] is the divergence of the tensor expression x,
which is the covariant derivative of x contracted on its last two
indices. If x is a k-tensor, Div[x] is a (k-1)-tensor. Div is the
formal adjoint of -Del.
\nOptions:
\n* Connection -> cn specifies that covariant derivatives are to be
taken with respect to the connection cn instead of the default
connection.
\n* Metric -> g: a symmetric 2-tensor expression, indicating that the
divergence is to be taken with respect to the Levi-Civita connection
of g instead of the default metric for x's bundle (which is assumed to
be Riemannian).";

DivGrad::usage = "DivGrad is a value for the global variable
$LaplacianConvention.";

Dot::usage = Ricci`Private`AppendMessage[Dot::usage,
"When x and y are tensor expressions, x.y represents the contraction of
the last index of x with the first index of y."];

Even::usage = "Even is a value for the Type option of DefineConstant.";

ERROR::usage = "If you attempt to plug in the wrong number of
indices in a tensor expression, Ricci will return ERROR[expression].";

Expand::usage = Ricci`Private`AppendMessage[Expand::usage,
"If you apply the Mathematica function Expand to an expression that
includes tensors with indices, Ricci automatically converts it to
TensorExpand."];

Extd::usage = "Extd[x] is the exterior
derivative of x, which must be an alternating
covariant tensor expression. In output form, Extd[x] prints as d[x].";

ExtdStar::usage = "ExtdStar[x] is the adjoint of the operator Extd
applied to x, which must be an alternating covariant tensor
expression.
\nOption:
\n* Metric -> g: a symmetric 2-tensor expression,
representing a metric to be used in place of the default
metric for x's bundle (which is assumed to be Riemannian).";

FactorConstants::usage = "FactorConstants[x] applies the Mathematica
function Factor to the constant factor in each term of the tensor
expression x. FactorConstants[x,n] or FactorConstants[x,{n1,n2,...}]
applies FactorConstants only to term n or to terms n1,n2,... of x.";

FirstBianchiRule::usage = "FirstBianchiRule is a rule that attempts to
turn sums containing two Riemannian curvature tensors into a single
term, using the first Bianchi identity.";

FirstStructureRule::usage = "FirstStructureRule is a rule that implements
the first structure equation for exterior derivatives of basis covectors.";

FirstUp::usage = "FirstUp is a value for the RiemannConvention option
of DefineBundle.";

FlatConnection::usage = "FlatConnection is a DefineBundle option.
\"FlatConnection -> True\" means the default connection has zero
curvature. The default is False.
\nThe function FlatConnection[bundle] returns True or False.";

FormQ::usage =
"FormQ[x] returns True if x is a Covariant alternating tensor
expression, and False otherwise."

Grad::usage =
"Grad[x] is the gradient of the tensor expression x. It is the same
as Del[x], except that the last index is Contravariant instead of
Covariant. If x is a function (0-tensor), then Grad[x] is a vector
field.";

Hermitian::usage =
"Hermitian is a value for the Symmetry option of DefineTensor. It is
valid only for 2-tensors.
A Hermitian
tensor is actually a Real, Symmetric 2-tensor associated with a bundle and
its conjugate, with the additional property that h[L[i],L[j]] and
h[LB[i],LB[j]] are both zero if i and j are associated with the same
bundle.";

HodgeInner::usage =
"HodgeInner[x,y] represents the Hodge inner product of the alternating
tensors x and y. In output form, HodgeInner[x,y]
appears as <<x, y>>.";

HodgeNorm::usage = "HodgeNorm[x] is the norm of the alternating
tensor expression x, with respect to the Hodge inner product. It
is automatically converted by Ricci to Sqrt[HodgeInner[x,Conjugate[x]]].";

Im::usage = Ricci`Private`AppendMessage[Im::usage,
"The Ricci package modifies Im to handle tensor expressions. Im[x]
is converted to (x - Conjugate[x])/(2I)."];

Imaginary::usage = "Imaginary is a value for the Type option of
DefineConstant, DefineTensor, and DefineMathFunction.";

IndexOrderedQ::usage = "IndexOrderedQ[{indices}] returns True if
the indices are ordered correctly according to Ricci's index ordering
rules:  first by name, then by altitude (lower before upper), and
False otherwise.";

IndexQ::usage = "IndexQ[i] returns True if i is an index name, and
False otherwise.";

Inner::usage = Ricci`Private`AppendMessage[Inner::usage, "If x and y are
tensor expressions of the same rank, Ricci interprets Inner[x,y] as the
inner product of x and y. In output form, Inner[x,y] appears as <x, y>."];

InsertIndices::usage = "If x is a tensor expression of rank k, then x[
i1,...,ik ] or InsertIndices[ x, {i1,...,ik} ] causes the indices
i1,...,ik to be inserted into the k index slots in order. If x is a
scalar (rank 0) expression, then InsertIndices[ x, {} ] or x[]
converts it to a component expression by generating dummy indices as
necessary, while InsertIndices[ x, {i1,...,ik} ] is automatically
replaced by CovD[ x, {i1,...,ik} ].";

Int::usage = "If x and y are alternating tensor expressions with
Rank[x] <= Rank[y], Int[x,y] represents the generalized interior
product of x into y. When x is a vector field, Int[x,y] is interior
multiplication of x into y (with a numerical factor depending on
$WedgeConvention). In general, Int[x,_] is the adjoint (with respect
to the Hodge inner product) of wedging with x on the left.";

Integer::usage = Ricci`Private`AppendMessage[Integer::usage, "Ricci
recognizes Integer as a value for the Type option of
DefineConstant."];

Inverse::usage = Ricci`Private`AppendMessage[Inverse::usage, "If x is a
2-tensor expression, Ricci interprets Inverse[x] as the inverse of x."];

Kronecker::usage =
"Kronecker[L[i],U[j]] is the Kronecker delta tensor.";

L::usage = "L[i] represents a lower index i.";

LaplaceBeltrami::usage = "LaplaceBeltrami[x] is the
Laplace-Beltrami operator applied to the differential form x.
It is automatically replaced by Extd[ExtdStar[x]] + ExtdStar[Extd[x]].
\nOption:
\n* Metric -> g: a symmetric 2-tensor expression, representing a
metric to be used in place of the default metric for x's bundle
(which is assumed to be Riemannian).";

Laplacian::usage = "Laplacian[x] is the covariant Laplacian of the
tensor expression x. It is automatically replaced by Div[Grad[x]] if
$LaplacianConvention = DivGrad (the default), and by -Div[Grad[x]] if
$LaplacianConvention = PositiveSpectrum.
\nOptions:
\n* Metric -> g: a symmetric 2-tensor expression, representing a
metric to be used in place of the default metric for the underlying
tangent bundle of x.
\n* Connection -> cn specifies that covariant derivatives are to be
taken with respect to the connection cn instead of the default
connection.";

LB::usage = "LB[i] represents a lower barred index i.";

LeviCivitaConnection::usage = "LeviCivitaConnection[g] represents the
Levi-Civita connection (as a tensor of rank {2,1}, i.e. a matrix of
one-forms) of the arbitrary metric g. When indices are inserted, the
components of the connection are computed in terms of the background
connection of g's bundle (assumed to be Riemannian) and covariant
derivatives of g.";

Lie::usage = "Lie[v,x] is the Lie derivative of the tensor expression
x in the direction v; v must be a vector field (a contravariant 1-tensor
expression).";

LieRule::usage = "LieRule transforms Lie derivatives of differential
forms to expressions involving exterior derivatives and Int.";

LowerAllIndices::usage = "LowerAllIndices[x] lowers all of the
indices in x by inserting appropriate metrics with raised indices.
LowerAllIndices[x,n] or LowerAllIndices[x,{n1,n2,...}] applies
LowerAllIndices only to term n or to terms n1,n2,... of x.";

Method::usage = Ricci`Private`AppendMessage[ Method::usage, "Method is an
option for the Ricci simplification command OrderDummy, which specifies
how hard the command should work to simplify the expression. The
allowable values are 0, 1, and 2."];

Metric::usage = "Metric[bundle] returns the bundle\'s metric.
\nMetric -> metric is an option for Del, Div, Grad, ExtdStar, CovD,
Laplacian, and LaplaceBeltrami.";

MetricQ::usage =
"MetricQ[x] returns True if x is a metric with or without indices,
and False otherwise.";

MetricType::usage = "MetricType is a DefineBundle option,
which specifies whether the bundle's metric has special properties.
Allowable values are MetricType -> Normal (no special
properties), and MetricType ->
Riemannian (for a Riemannian tangent bundle). The
default is Normal. If Riemannian is specified, the following additional
options may also be given:
\n* RiemannTensor -> name: a name to
be given to the Riemannian curvature tensor
for this bundle. Default is Rm.
\n* RicciTensor -> name: a name to be given to the Ricci tensor
for this bundle. Default is Rc.
\n* ScalarCurv -> name: a name to be given to the scalar curvature function
for this bundle. Default is Sc.
\n* RiemannConvention -> SecondUp or FirstUp: determines the sign convention
used for the Riemannian curvature tensor for this bundle. Default is SecondUp,
or the value of the global variable $RiemannConvention if it has been set.";

Negative::usage = Ricci`Private`AppendMessage[Negative::usage, "Ricci
recognizes Negative as a value for the Type option of DefineConstant,
DefineTensor, DefineMathFunction."];

NewDummy::usage = "NewDummy[x] converts all dummy indices occurring in
the tensor expression x to computer-generated dummy indices.
NewDummy[x,n] or NewDummy[x,{n1,n2,...}] applies NewDummy only to term
n or to terms n1,n2,... of x. For most bundles, the dummy names are
of the form \"kn\", where k=BundleDummy[bundle] and n is an integer.
For one-dimensional bundles, only the dummy name k itself is
generated.";

NewRule::usage = "NewRule is a DefineRule option. NewRule -> True
specifies that previous rules defined with the same name are to be
erased before defining the new rule. The default is False, which
means that this rule is to be appended to the list of already-existing
rules by the same name.";

NonNegative::usage = Ricci`Private`AppendMessage[NonNegative::usage, "Ricci
recognizes NonNegative as a value for the Type option of DefineConstant,
DefineTensor, DefineMathFunction."];

NonPositive::usage = Ricci`Private`AppendMessage[NonPositive::usage, 
"NonPositive is a value for the Type option of
DefineConstant, DefineTensor, DefineMathFunction."];

NoOneDims::usage = "NoOneDims is a value for the Mode option of
AbsorbMetrics.";

Norm::usage =  Ricci`Private`AppendMessage[Norm::usage,"Norm[x]
is the norm of the tensor expression x. It
is automatically converted by Ricci to Sqrt[Inner[x,Conjugate[x]]]."];

Normal::usage = Ricci`Private`AppendMessage[Normal::usage, 
"Normal is a value for the MetricType option of
DefineBundle."];

NoSymmetries::usage = "NoSymmetries is a value for the
Symmetries option of DefineTensor.";

Odd::usage = "Odd is a value for the Type option of DefineConstant.";

OneDims::usage = "OneDims is an AbsorbMetrics option. OneDims -> False
means that AbsorbMetrics should not absorb one-dimensional metrics
(unless they are contracted with other metrics). Default is True.";

OrderCovD::usage = "OrderCovD[x] orders all of the indices appearing
after \";\" in the tensor expression x, by adding appropriate
curvature and torsion terms. OrderCovD[x,n] or
OrderCovD[x,{n1,n2,...}] applies OrderCovD only to term n or to terms
n1,n2,... of x. ";

OrderDummy::usage = "OrderDummy[x] attempts to put the dummy indices
occurring in the tensor expression x in a \"canonical form\".
OrderDummy[x,n] or OrderDummy[x,{n1,n2,...}] applies OrderDummy only
to term n or to terms n1,n2,... of x. Dummy index pairs are ordered
so that the lower member appears first whenever possible. OrderDummy
tries various rearrangements of dummy index names, and chooses the
lexically smallest version of the expression that results." <>
"\nOption:
\n* Method -> n specifies how hard OrderDummy should work to find the
best possible version of the expression. The default is Method -> 1,
which means that dummy index names will be interchanged in pairs only.
Method -> 2 causes OrderDummy to try all possible permutations of the
dummy index names. Method -> 0 means don't try interchanging names at
all.";

OrthonormalFrame::usage = "OrthonormalFrame is an option for
DefineBundle. If True, it specifies that the metric coefficients are
always assumed to be constants in the default basis. For
one-dimensional bundles, the metric coefficient is always taken to be
1. Default is False";

ParallelFrame::usage = "ParallelFrame is an option for DefineBundle.
If True, it specifies that the default basis is always assumed to be
parallel. Thus the connection and curvature forms will always be
zero, and the bundle will be flat. Default is False.";

Plus::usage = Ricci`Private`AppendMessage[Plus::usage, "Ricci
transforms expressions of the form (a + b)[L[k],...] into
InsertIndices[a+b,{L[k],...}]."];

Positive::usage = Ricci`Private`AppendMessage[Positive::usage, "Ricci
recognizes Positive as a value for the Type option of DefineConstant,
DefineTensor, DefineMathFunction."];

PositiveDefinite::usage = "PositiveDefinite is a DefineBundle option.
PositiveDefinite -> False means that the bundle's metric is not assumed
to be positive definite. Default is True.";

PositiveInteger::usage = "PositiveInteger is an obsolete value for the
Type option of DefineConstant. Use Type -> {Positive,Integer}
instead.";

PositiveInteger := (Message[PositiveInteger::usage];{Positive,Integer});

PositiveReal::usage = "PositiveReal is an obsolete value for the Type
option of DefineConstant, DefineTensor, and DefineMathFunction. Use
Type -> {Positive,Real} instead.";

PositiveReal := (Message[PositiveReal::usage];{Positive,Real});

PositiveSpectrum::usage = "PositiveSpectrum is a value for the global
variable $LaplacianConvention.";

Power::usage = Ricci`Private`AppendMessage[Power::usage, "If x is a
tensor expression and p is a positive integer, Ricci interprets x^p as
the p-th symmetric power of x with itself. \nRicci transforms
expressions of the form (x^p)[L[i],...] into
InsertIndices[x^p,{L[i],...}] or CovD[x^p,{L[i],...}]." <> "\nRicci
causes a power of a product such as (a b)^p to be expanded into a
product of powers a^p * b^p, provided a and b do not contain any
indices; if they do contain indices, then (a b)^p is transformed to
Summation[a b]^p, to prevent the expression from being expanded to a^p
b^p. Summation is not printed in output form."];

PowerSimplify::usage = "PowerSimplify[x] attempts to simplify negative and
nonintegral powers that appear in the tensor expression x, by
expanding and collecting constants in the base and the exponent of
each such power. PowerSimplify[x,n] or PowerSimplify[x,{n1,n2,...}] applies
PowerSimplify only to term n or to terms n1,n2,... of x.";

ProductExpand::usage = "ProductExpand[x] expands out symmetric
products and wedge products of 1-tensors that occur in x, and rewrites
them in terms of tensor products.";

Quiet::usage = "Quiet is an option for some of the defining and
undefining commands in the Ricci package. The option Quiet -> True will
silence the usual messages printed by these commands. The default is
the value of the global variable $Quiet, which is initially False."

Rank::usage = "Rank[x] returns the rank of the tensor expression x.";

Re::usage = Ricci`Private`AppendMessage[Re::usage,
"The Ricci package modifies Re to handle tensor expressions.
Re[x] is converted to (x + Conjugate[x])/2."];

Real::usage = Ricci`Private`AppendMessage[Real::usage, "Ricci
recognizes Real as a value for the Type option of DefineBundle,
DefineConstant, DefineTensor, and DefineMathFunction."];

RenameDummy::usage = "RenameDummy[x] changes the names of dummy
indices in x to standard names. RenameDummy[x,n] or
RenameDummy[x,{n1,n2,...}] applies RenameDummy only to term n or to
terms n1,n2,... of x. RenameDummy chooses names in alphabetical order
from the list of index names associated with the appropriate bundle,
skipping those names that already appear in x as free indices. When
the list of index names is exhausted, computer-generated names of the
form \"kn\" are used, where k is the last index name in the list and n
is an integer. For one-dimensional bundles, n is omitted.";

Ricci::usage = "Ricci is a Mathematica package for doing symbolic
computations that arise in differential geometry. To load it, type
\"<<Ricci.m\". To find out more information about a Ricci command,
object, or function, type \"?name\".";

RicciSave::usage = "RicciSave[\"filename\"] causes all of your current Ricci
definitions to be saved into the file \"filename\". The previous
contents of \"filename\" are erased. The definitions can be read back
in by the command <<filename.";

RicciTensor::usage = "RicciTensor -> name is a DefineBundle option,
which specifies a name to be used for the Ricci tensor for this
bundle. The default is Rc. This option takes effect only when
MetricType -> Riemannian has been specified.
\nRicciTensor[g] represents the Ricci curvature of the arbitrary
metric g.";

RiemannConvention::usage = "RiemannConvention -> SecondUp or FirstUp
is a DefineBundle option, which specifies the index convention of the
Riemannian curvature tensor. The default is SecondUp, or the value of
the global variable $RiemannConvention if it has been set. This
option takes effect only when MetricType -> Riemannian has been
specified.";

Riemannian::usage =
"MetricType -> Riemannian is a DefineBundle option for defining a
Riemannian tangent bundle. If this option
is specified, the following additional
options may be given:
\n* RiemannTensor -> name: the name for the Riemannian curvature
tensor for this bundle. Default is Rm.
\n* RicciTensor -> name: the name for the Ricci curvature
for this bundle. Default is Rc.
\n* ScalarCurv -> name: the name for the scalar
curvature for this bundle. Default is Sc.
\n* RiemannConvention -> FirstUp or SecondUp: controls the index convention
of the Riemannian curvature tensor. Default is SecondUp.";

RiemannSymmetries::usage = "RiemannSymmetries is a value for the
Symmetries option of DefineTensor. Tensors having RiemannSymmetries
must be of degree 4. The symmetries are those of the Riemann
curvature tensor, namely skew-symmetric in the first two indices,
skew-symmetric in the last two indices, symmetric when the first two
and last two are interchanged. The (first) Bianchi identity can be
explicitly applied to any tensor with these symmetries by means of
FirstBianchiRule or BianchiRules.";

RiemannTensor::usage = "RiemannTensor -> name is a DefineBundle
option, which specifies a name to be used for the Riemannian curvature
tensor for this bundle. The default is Rm. This option takes effect
only when MetricType -> Riemannian has been specified.
\nRiemannTensor[g] represents the Riemannian curvature tensor of the
arbitrary metric g.";

Same::usage = "Bundle -> Same is a special DefineTensor option for
product tensors. It means that indices from any bundle can be
inserted in the first set of index slots, but indices in the remaining
slots must be from the same bundle or its conjugate.";

ScalarCurv::usage = "ScalarCurv -> name is a DefineBundle option,
which specifies a name to be used for the scalar curvature for this
bundle. The default is Sc. This option takes effect only when
MetricType -> Riemannian has been specified.
\nScalarCurv[g] represents the scalar curvature of the arbitrary
metric g.";

ScalarQ::usage = "ScalarQ[x] is True if x is a rank 0 tensor
expression with no free indices, and False otherwise.";

SecondBianchiRule::usage = "SecondBianchiRule is a rule that attempts
to turn sums containing two differentiated Riemannian curvature
tensors into a single term using the second Bianchi identity.";

SecondStructureRule::usage = "SecondStructureRule is a rule that
implements the second structure equation for exterior derivatives of
the generic connection forms.";

SecondUp::usage = "SecondUp is a value for the RiemannConvention option
of DefineBundle.";

SimplifyConstants::usage = "SimplifyConstants[x] applies the
Mathematica function Simplify to the constant factor in each term of
the tensor expression x. SimplifyConstants[x,n] or
SimplifyConstants[x,{n1,n2,...}] applies SimplifyConstants only to
term n or to terms n1,n2,... of x.";

Skew::usage = "Skew is a synonym for Alternating.";

SkewHermitian::usage = "SkewHermitian is a value for the Symmetry
option of DefineTensor. A SkewHermitian tensor is actually a Real,
Alternating tensor associated with a bundle and its conjugate, with
the additional property that h[L[i],L[j]] and h[LB[i],LB[j]] are both
zero if i and j are associated with the same bundle.";

SkewQ::usage = "SkewQ[x] is True if x is an alternating or
skew-Hermitian tensor expression, and False otherwise.";

StructureRules::usage = "StructureRules is the union of
CompatibilityRule, FirstStructureRule, and SecondStructureRule.";

Summation::usage = "If a and b contain indices, then Ricci transforms (a
b)^p internally to Summation[a b]^p, to prevent the expression from
being expanded to a^p b^p. Summation is not printed in output
form.";

SuperSimplify::usage = "SuperSimplify[x] attempts to put the tensor
expression x into a canonical form, so that two expressions that are
equal will usually be identical. SuperSimplify[x,n] or
SuperSimplify[x,{n1,n2,...}] applies SuperSimplify only to term n or
to terms n1,n2,... of x. SuperSimplify works exactly the same way as
TensorSimplify, except that it calls OrderDummy with the option Method
-> 2, so that all possible permutations of dummy index names will be
tried. This will generally be much slower for expressions having more
than 4 dummy index pairs.";

Sym::usage = "Sym[x] represents the symmetrization of the tensor
expression x.";

Symmetric::usage = "Symmetric is a value for the Symmetries option
of DefineTensor.";

SymmetricProduct::usage = "If x, y, and z are tensor expressions,
SymmetricProduct[x,y,z] or x * y * z or x y z represents their symmetric
product. Symmetric products are represented internally by ordinary
multiplication. In output form, Ricci inserts explicit asterisks
whenever non-scalar tensor expressions are multiplied together, to
remind the user that multiplication is being interpreted as symmetric
product.";

SymmetricQ::usage = "SymmetricQ[x] is True if x is a symmetric or
Hermitian tensor expression, and False otherwise.";

Symmetries::usage = "Symmetries -> sym is a DefineTensor option.
Allowed values are NoSymmetries, Symmetric, Alternating, Hermitian,
SkewHermitian, RiemannSymmetries, or a symmetry defined with
DefineTensorSymmetries. Skew is a synonym for Alternating. Hermitian
and SkewHermitian can be used only with 2-tensors. For product
tensors, sym can be a list of symmetries, one for each rank in the
rank list.";

TangentBundle::usage = "TangentBundle is a DefineBundle option. If
bundle x is defined with the option \"TangentBundle -> {y,z}\", this
means that the underlying tangent bundle of x is the direct sum of y
and z. The form \"TangentBundle -> y\" may be used when there is only
one bundle. If $DefaultTangentBundle is defined, it will be used as
the default when TangentBundle is not specified; otherwise by default
the tangent bundle of x is assumed to be x itself.
\nThe function
TangentBundle[x] returns the tangent bundle list for the bundle x.";

Tensor::usage = "Internally, Ricci represents tensors in the form
Tensor[name,{i,j,...},{k,l,...}], where i,j,... are the tensor indices
and k,l,... are indices resulting from covariant differentiation. In
input form, an unindexed tensor is represented just by typing its
name. Indices are inserted by typing them in brackets after the
tensor name. Once all index slots are full, indices resulting from
covariant differentiation can be typed in a second set of brackets.
Tensors are created by DefineTensor and removed by UndefineTensor.";

TensorCancel::usage = "TensorCancel[x] attempts to simplify each term
of x by canceling common factors, even when the factors have
different names for their dummy indices. TensorCancel[x,n] or
TensorCancel[x,{n1,n2,...}] applies TensorCancel only to term n or to
terms n1,n2,... of x.";

TensorData::usage = "TensorData[name] is a list containing data for
the tensor \"name\", used internally by Ricci.";

TensorExpand::usage = "TensorExpand[x] expands products and
positive integral powers in x, just as Expand does, but maintains
correct dummy index conventions and does not expand constant factors.
TensorExpand[x,n] or TensorExpand[x,{n1,n2,...}] applies TensorExpand
only to term n or to terms n1,n2,... of x.";

TensorFactor::usage = "TensorFactor[expr] returns the product of all the
non-constant factors in expr, which should be a monomial.";

TensorMetricQ::usage = "TensorMetricQ[tensorname] returns True if
tensorname is the metric of some bundle, and False otherwise.";

TensorProduct::usage =
"TensorProduct[x,y,z] or TProd[x,y,z] or x ~TProd~ y ~TProd~ z
represents the tensor product of x, y, and z.
In output form, tensor products appear as
\n   x (X) y (X) z.";

TensorQ::usage = "TensorQ[name] returns True if name is the name of a
tensor, and False otherwise.";

TensorRankList::usage = "TensorRankList stores the list of ranks for a
product tensor. Used internally by Ricci.";

TensorSimplify::usage = "TensorSimplify[x] attempts to put the tensor
expression x into a canonical form, so that two expressions that are
equal will usually be identical. TensorSimplify[x,n] or
TensorSimplify[x,{n1,n2,...}] applies TensorSimplify only to term n or
to terms n1,n2,... of x. TensorSimplify expands products and positive
integer powers, uses metrics to raise and lower indices, tries to
rename all dummy indices in a canonical order, and collects all terms
containing the same tensor factors but different constant factors. It
does not reorder indices after the \";\" (use OrderCovD or
CovDSimplify to do that). TensorSimplify calls CorrectVariances,
TensorExpand, AbsorbMetrics, PowerSimplify, RenameDummy, OrderDummy, and
CollectConstants.";

TensorSymmetry::usage = "A tensor symmetry is an object of the form
TensorSymmetry[name,d,{perm1,sgn1,...,permN,sgnN}], where d is a positive
integer, permi is a non-trivial permutation of {1,2,...,d}, and sgni is
plus or minus 1. TensorSymmetry objects can be defined with
DefineTensorSymmetries and used in the
DefineTensor command.";

TeXFormat::usage = "TeXFormat -> \"string\" is an option for DefineTensor
and DefineIndex. It specifies the way the tensor or index name will
appear in TeXForm. The default is the tensor's or index's name.";

Times::usage = Ricci`Private`AppendMessage[Times::usage,
"Ricci uses ordinary multiplication
to represent symmetric products of tensors.
Ricci transforms an
expression of the form (a * b)[i]
into InsertIndices[a*b,{i}] whenever \"i\" is
one or more indices (L[j] or U[j]).
\nIn output form, Ricci modifies Mathematica's usual
ordering of factors: constants are printed first, then fully-indexed tensors,
then other tensor expressions."];

Tor::usage = "Tor is the name for the generic torsion tensor for the
default connection in any bundle. It is a product tensor of rank
{2,1}; inserting the first two indices yields the torsion 1-forms
associated with the default basis. Inserting all three indices yields
the components of the torsion tensor.";

TorsionFree::usage = "TorsionFree is a DefineBundle option.
\"TorsionFree -> False\" means the default connection has nonvanishing
torsion. The default is True.
\nThe function TorsionFree[bundle]
returns True or False.";

TotalRank::usage = "TotalRank[tensorname] returns the total rank of
the tensor tensorname.";

TProd::usage = "TProd is an abbreviation for TensorProduct.";

Tr::usage = Ricci`Private`AppendMessage[ Tr::usage,
"If x is a 2-tensor expression, Ricci interprets
Tr[x] as the trace of the x with respect to the metric of x's bundle." ];

Transpose::usage = Ricci`Private`AppendMessage[ Transpose::usage, "If x
is a tensor expression, Ricci interprets Transpose[x] as x with its
index positions reversed." ];

Type::usage = "Type is an option for DefineBundle, DefineConstant,
DefineTensor, and DefineMathFunction.";

U::usage = "U[i] represents an upper index i.";

UB::usage = "UB[i] represents an upper barred index i.";

UndefineBundle::usage = "UndefineBundle[bundle] clears the definition
of bundle. It also clears the definition of the bundle's metric and
indices.
\nOption:
\n* Quiet -> True or False. Default is False, or the value of $Quiet if set.";

UndefineConstant::usage =
"UndefineConstant[symbol] removes symbol's definition as a constant.
\nOption:
\n* Quiet -> True or False. Default is False, or the value of $Quiet if set.";

UndefineIndex::usage =
"UndefineIndex[i] or UndefineIndex[{i,j,k}] removes the association
of indices with their bundles.
\nOption:
\n* Quiet -> True or False. Default is False, or the value of $Quiet if set.";

UndefineRelation::usage =
"UndefineRelation[tensor] deletes the relation previously defined for tensor.
The tensor must exactly match the first argument of the corresponding call
to DefineRelation.
\nOptions:
\n* Quiet -> True or False. Default is False, or the value of $Quiet if set.";

UndefineTensor::usage = "UndefineTensor[tensorname] clears the definition
of tensorname.
\nOption:
\n* Quiet -> True or False. Default is False, or the value of $Quiet
if set.";

UndefineTensorSymmetries::usage = "UndefineTensorSymmetries[name] deletes
the TensorSymmetry object created by DefineTensorSymmetries[name,...].";

UnderlyingTangentBundle::usage = "UnderlyingTangentBundle[x] returns a
list of bundles representing the underlying tangent bundle of the
expression. The tangent bundle is the direct sum of the bundles in
the list. It is assumed that all tensors used in a given expression
have the same underlying tangent bundle.";

Variance::usage = Ricci`Private`AppendMessage[ Variance::usage,
"Variance -> Covariant or Contravariant is a
DefineTensor option. May be a list whose length is equal to the total
rank of the tensor, in which each entry specifies the variance of the
corresponding index slot. May be abbreviated Co and Con. Default is
Covariant.
\nIf x is a tensor expression, Variance[x] returns a list of variances of
x, one for each index slot." ];

VectorFieldQ::usage = "VectorFieldQ[x] returns True if x is a
Contravariant 1-tensor expression, and False otherwise.";

Wedge::usage = "If x, y, and z are alternating tensor expressions,
then x ~Wedge~ y ~Wedge~ z or Wedge[x,y,z] represents their wedge
product. The interpretation of Wedge in terms of tensor products is
determined by the global variable $WedgeConvention. In output form,
wedge products print as x ^ y ^ z.";

(************************* Global variables ***********************)

$DefaultTangentBundle::usage = "The global variable
$DefaultTangentBundle can be set by the user to a bundle or list of
bundles. It is used by DefineBundle as the default value for the
TangentBundle option, and by DefineTensor as the default value for the
Bundle option. By default, it is set to the first bundle the user
defines (or the direct sum of this bundle and its conjugate if the bundle
is complex).";

$LaplacianConvention::usage = "$LaplacianConvention determines which sign
convention is used for the covariant Laplacian on functions and tensors.
$LaplacianConvention = DivGrad means that Laplacian[x] = Div[Grad[x]],
while $LaplacianConvention = PositiveSpectrum means that Laplacian[x] =
-Div[Grad[x]]. Default is DivGrad.";

$MathFunctions::usage = "$MathFunctions is a list of names that have
been defined as scalar mathematical functions for use by Ricci.";

$Quiet::usage = "The global variable $Quiet is used by the Ricci package
to determine whether the defining and undefining commands report on
what they are doing. Setting $Quiet=True will silence these chatty
commands. The default is False. It can be overridden for a
particular command by specifying the Quiet option as part of the
command call.";

$RiemannConvention::usage = "The global variable $RiemannConvention
can be set by the user to specify a default value for the
RiemannConvention option of DefineBundle and MakeRiemannian. Default
is SecondUp.";

$TensorFormatting::usage = "The global variable $TensorFormatting can
be set to True or False by the user to turn on or off Ricci's special
output formatting of tensors and indices. Default is True.";

$TensorTeXFormatting::usage = "The global variable
$TensorTeXFormatting can be set to True or False by the user to turn
on or off Ricci's special formatting of tensors in TeXForm. Default is
True.";

$WedgeConvention::usage = "The global variable $WedgeConvention can be
set by the user to determine the interpretation of wedge products.
The allowed values are Alt and Det. The default is Alt.";

(********************* Error messages *************************)

General::invalid = "\"`1`\" is not a valid `2`.";

Index::error = "Incorrect number of indices provided for `1`.";

Bundle::error = "No bundle specified or implied for `1`.";

General::oneterm = "`1` has only one term.";

General::termspec = "`1` is not a valid term number or list
of term numbers.";

(************ Initialize default values for user parameters *************)

$Quiet=False;
$DefaultTangentBundle = Null;
$WedgeConvention = Alt;
$RiemannConvention = SecondUp;
$TensorFormatting = True;
$TensorTeXFormatting = True;
$MathFunctions = {};
$LaplacianConvention = DivGrad;

(******************** Switch to private context ************************)

Begin["`Private`"];

Unprotect[TensorQ];
SetAttributes[ { TensorQ,Tensor,TensorData,TensorMetricQ,TotalRank,rank }, 
 HoldFirst];
SetAttributes[ ERROR, HoldAll ];
Protect[TensorQ];

(********************** General utility routines **************************)

(** MakeLinear can be used to make a function behave linearly with respect
    to sums and constant multiples. **)

MakeLinear[f_Symbol] := (
  HoldPattern[ f[(c_?ConstantQ) * v_] ] := c f[v];
  HoldPattern[ f[x_Plus] ] := f /@ x;
  );

(* PartAt[ exp, pos ] behaves like Part[exp,pos], except that "pos" is
   a list.  *)

PartAt[ exp_, {pos___Integer} ] := Part[exp,pos];
HeldPartAt[ exp_, {pos___Integer} ] := HeldPart[exp,pos];

(* Occurrences[exp,pattern] returns a list of the actual occurrences
   of "pattern" in exp. *)

Occurrences[ exp_, pattern_ ] := PartAt[ exp, # ]& /@ Position[exp,pattern];

(* ContainsQ[list1,list2] is True if list2 is a subset of list1, 
   False otherwise. *)

ContainsQ[ a_, b_ ] := And @@ (MemberQ[a,#]& /@ b);

(** optionsOK is an internal subroutine to check the validity of optional
   arguments.  optionsOK[command, options] returns True if the "options"
   represents a valid list of options for "command", and prints a message
   and returns False if not. **)

optionsOK[subr_,opts___] := Module[{badoptions},
 badoptions = Select[{opts},
   !(Head[#]===Rule && MemberQ[optionNames[Options[subr]],#[[1]]])&];
 If[badoptions==={},
   (*then OK*)    Return[True],
   (*else error*) Message[
                    Evaluate[ToExpression[ToString[subr]<>"::invalid"]],
                    badoptions[[1]],
                    ToString[subr]<>" option"];
                  Return[False]
   ]];

(** optionNames takes a list of rules and gives back the list of
    left-hand sides.  For example:
    optionNames[{a->b,c->d}] = {a,c}.
**)

optionNames[optionlist_] := #[[1]]& /@ optionlist;

(* newsymQ checks that a symbol hasn't been used before as a tensor or bundle
   name, and that it hasn't been assigned a value. *)

SetAttributes[newsymQ,HoldFirst];

newsymQ[Conjugate[x_]] := newsymQ[x];
newsymQ[x_] :=
  (!ValueQ[x] && Head[x]===Symbol && Attributes[x]==={} &&
   !BundleQ[x] && !IndexQ[x]);

tfQ[True|False] := True;
tfQ[_] := False;

(************************** RicciSave **********************************)

RicciSave[] := Message[RicciSave::argx,"RicciSave",0];
RicciSave[_,extras__] := Message[RicciSave::argx,"RicciSave",Length[{extras}]+1];

RicciSave[file_Symbol] := RicciSave[ToString[file]];

RicciSave[file_String] := Module[{tform,protectlist},
  Put @@ Append[ deflist["Ricci`$"], file ];
  tform = $TensorFormatting;
  $TensorFormatting = False;
  protectlist = Unprotect[Evaluate[namelist[$Context]]];
  PutAppend @@ Append[ deflist[$Context], file ];
  PutAppend[ "DefineMathFunction[$MathFunctions,Quiet->True]"//OutputForm,
             file ];
  Unprotect[ Curv,Tor,Conn,Basis,Kronecker];
  PutAppend[ "Unprotect[Curv,Tor,Conn,Basis,Kronecker]"//OutputForm, file];
  PutAppend[ Definition[Curv,Tor,Conn,Basis,Kronecker], file ];
  PutAppend[ "Protect[Curv,Tor,Conn,Basis,Kronecker]"//OutputForm, file];
  PutAppend[ "Protect["//OutputForm, file];
  PutAppend[ protectlist, file];
  PutAppend[ "]"//OutputForm, file];
  $TensorFormatting = tform;
  Protect[ Curv,Tor,Conn,Basis,Kronecker];
  Protect[ Evaluate[protectlist]];
  ];

RicciSave[x_] := Message[RicciSave::invalid, x, "file name"];

namelist[prefix_] := Names[ StringJoin[prefix, "*"] ];

deflist[prefix_] :=
  Select[ Definition /@ namelist[prefix],
          ToString[#] =!= "Null" &];

(*************************** Bundle *****************************)

SetAttributes[{TensorMetricQ,DefineBundle},HoldFirst];

(* MetricQ[exp] gives True iff exp is a metric 2-tensor with or without
   indices.  If exp is a tensor, this is translated into 
   TensorMetricQ, which is a HoldFirst function so that it
   can be associated with the name of the tensor.  *)

MetricQ[Tensor[g_,{_,_},{}]] := TensorMetricQ[g];
MetricQ[Tensor[g_,{},{}]] := TensorMetricQ[g];
MetricQ[_] = False;

(* Initialize TensorMetricQ[anything] to be False. 
   DefineBundle will set this to True for any tensors that are metrics.  *)

TensorMetricQ[_] = False;

(* BundleQ[name] = True iff name is a bundle.  Initially false for any 
   name, it is set to True by DefineBundle. *)

BundleQ[_] = False;

(* Tell the various bundle information functions how to deal with conjugate
   bundles. *)
BundleQ[Conjugate[b_]] := BundleQ[b];
Dimension[Conjugate[b_]] := Dimension[b];
Metric[Conjugate[b_]] := Metric[b];
TorsionFree[Conjugate[b_]] := TorsionFree[b];
FlatConnection[Conjugate[b_]] := FlatConnection[b];
ParallelFrame[Conjugate[b_]] := ParallelFrame[b];
OrthonormalFrame[Conjugate[b_]] := OrthonormalFrame[b];
CommutingFrame[Conjugate[b_]] := CommutingFrame[b];
Type[Conjugate[b_]] := Type[b];
BundleIndices[Conjugate[b_]] := Conjugate /@ BundleIndices[b];
TangentBundle[Conjugate[b_]] := TangentBundle[b];
PositiveDefinite[Conjugate[b_]] := PositiveDefinite[b];

Dimension[Any] := Plus @@ Dimension /@ Flatten[{$DefaultTangentBundle}];

(* These are needed by UnderlyingTangentBundle. *)
TangentBundle[Any] := $DefaultTangentBundle;
TangentBundle[Null] := {Null};

(* The following are needed when computing conjugates of bundles *)

Unprotect[{Null,Automatic}];
Any/:       Conjugate[Any] := Any;
Null/:      Conjugate[Null] := Null;
Automatic/: Conjugate[Automatic] := Automatic;
Protect[{Null,Automatic}];

(* Initialize BundleIndices[any bundle] to an empty list *)
BundleIndices[_] := {};

(* BundleDummy returns the last index in BundleIndices *)
BundleDummy[bun_] := Last @ BundleIndices @ bun;

(*********** argument-checking subroutines for DefineBundle **********)

dimQ[d_] := ConstantQ[d] && (!NumberQ[d] || (IntegerQ[d] && d > 0));

decompQ[Automatic] = True;
decompQ[Null] = True;
decompQ[x_List] := And @@ ((BundleQ[#] || newsymQ[#])& /@ x);
decompQ[x_] := decompQ[{x}];

indQ[{i__},name_]:= And @@ (indQ[#,name]&) /@ {i}; 
indQ[i_Symbol,name_] := newsymQ[i] || (IndexQ[i] && Bundle[i]===name);
indQ[_,_] := False;

mtypeQ[t_] := t===Normal || t===Riemannian;

rmconvQ[FirstUp] = True;
rmconvQ[SecondUp] = True;
rmconvQ[_] = False;

(* Define the default names for Riemannian curvature tensors.
   These should not be evaluated until used, so that the names
   are defined in the user's context. *)

defaultRm := ToExpression["Rm"];
defaultRc := ToExpression["Rc"];
defaultSc := ToExpression["Sc"];

(***************** DefineBundle ***********************************)

Options[DefineBundle] :=
  {Type -> Real,
   TangentBundle -> Automatic,
   FlatConnection -> False,
   TorsionFree -> True,
   Quiet -> $Quiet,
   RiemannTensor -> Null,
   RicciTensor -> Null,
   ScalarCurv -> Null,
   RiemannConvention -> $RiemannConvention,
   MetricType -> Normal,
   BasisName -> Automatic,
   CoBasisName -> Automatic,
   ParallelFrame -> False,
   OrthonormalFrame -> False,
   CommutingFrame -> False,
   PositiveDefinite -> True};

DefineBundle[]      := Message[DefineBundle::argm,"DefineBundle",0,4];
DefineBundle[_]     := Message[DefineBundle::argmu,"DefineBundle",4];
DefineBundle[_,_]   := Message[DefineBundle::argm,"DefineBundle",2,4];
DefineBundle[_,_,_] := Message[DefineBundle::argm,"DefineBundle",3,4];

DefineBundle[name_, dim_, metric_, indices_, opts___] := 
  Module[{type,compat,tdecomp,flat,tfree,quiet,rm,rc,sc,rmc,mtype,
          riemanntensor,riccitensor,scalarcurv,bname,cobname,pframe,
          oframe,cframe,pdef},
  If [!optionsOK[DefineBundle,opts], Return[]];
  {type,tdecomp,flat,tfree,quiet,rm,rc,sc,
   rmc,mtype,bname,cobname,pframe,oframe,cframe,pdef} =
    {Type,TangentBundle,FlatConnection,TorsionFree,Quiet,
     RiemannTensor,RicciTensor,ScalarCurv,RiemannConvention,
     MetricType,BasisName,CoBasisName,ParallelFrame,
     OrthonormalFrame,CommutingFrame,PositiveDefinite} /.
    {opts} /. Options[DefineBundle];
  Which[
   !newsymQ[name], 
     Message[DefineBundle::invalid,HoldForm[name],"bundle name"],
   !newsymQ[metric],
     Message[DefineBundle::invalid,metric,"metric name"],
   !dimQ[dim],
     Message[DefineBundle::invalid,dim,"dimension"],
   !indQ[indices,name],
     Message[DefineBundle::invalid,indices,"list of indices"],
   dim===1 && !Length[Flatten[{indices}]]===1,
     Print["ERROR: a one-dimensional bundle can have only one index name"],
   !MemberQ[{Real,Complex}, type],
     Message[DefineBundle::invalid,type,"Type"],
   !decompQ[tdecomp],
     Message[DefineBundle::invalid,tdecomp,"tangent bundle"],
   !tfQ[flat],
     Message[DefineBundle::invalid,flat,"value for FlatConnection"],
   !tfQ[tfree],
     Message[DefineBundle::invalid,tfree,"value for TorsionFree"],
   !mtypeQ[mtype],
     Message[DefineBundle::invalid,mtype,"MetricType"],
   !rmconvQ[rmc],
     Message[DefineBundle::invalid,rmc,"RiemannConvention"],
   rm =!= Null && !newsymQ[Evaluate[rm]],
     Message[DefineBundle::invalid,rm,"RiemannTensor name"],
   mtype === Riemannian && rm === Null && !newsymQ[Evaluate[defaultRm]],
     Message[DefineBundle::invalid,defaultRm,"RiemannTensor name"],
   rc =!= Null && !newsymQ[Evaluate[rc]],
     Message[DefineBundle::invalid,rc,"RicciTensor name"],
   mtype === Riemannian && rc === Null && !newsymQ[Evaluate[defaultRc]],
     Message[DefineBundle::invalid,defaultRc,"RicciTensor name"],
   sc =!= Null && !newsymQ[Evaluate[sc]],
     Message[DefineBundle::invalid,sc,"ScalarCurv name"],
   mtype === Riemannian && sc === Null && !newsymQ[Evaluate[defaultSc]],
     Message[DefineBundle::invalid,defaultSc,"ScalarCurv name"],
   bname =!= Automatic && !newsymQ[Evaluate[bname]],
     Message[DefineBundle::invalid,bname,"BasisName"],
   cobname =!= Automatic && !newsymQ[Evaluate[bname]],
     Message[DefineBundle::invalid,cobname,"CoBasisName"],
   !tfQ[pframe],
     Message[DefineBundle::invalid,pframe,"value for ParallelFrame"],
   !tfQ[oframe],
     Message[DefineBundle::invalid,oframe,"value for OrthonormalFrame"],
   !tfQ[cframe],
     Message[DefineBundle::invalid,cframe,"value for CommutingFrame"],
   !tfQ[pdef],
     Message[DefineBundle::invalid,pdef,"value for PositiveDefinite"],
   (* Now all the arguments are OK.  Let's define the bundle. *)
   True,
     name/: BundleQ[name] = True;
     If[ type===Real, (*then*) (name/: Conjugate[name] = name)];
     If[ newsymQ[dim],
         DefineConstant[dim,Type->{Positive,Integer}] ];
     name/: Dimension[name] = dim;
     (* If $DefaultTangentBundle hasn't been set, then set it as best
        we can from the information available. *)
     Which[ $DefaultTangentBundle =!= Null,
              Null, (* do nothing *)
            tdecomp =!= Automatic,
              $DefaultTangentBundle = tdecomp,
            True,
              $DefaultTangentBundle = Union[{name,Conjugate[name]}] ];
     (* Flatten appears below to make sure we get a list, even if
        the user just specified a single bundle *)
     name/: TangentBundle[name] = Flatten[{tdecomp} /.
                                  Automatic -> $DefaultTangentBundle ];
     name/: Type[name] = type;
     (* Associate the indices with the bundle *)
     DefineIndex[indices,name,Quiet->quiet];
     (* Define the metric as a symmetric or hermitian 2-tensor *)
     DefineTensor[metric,2,Type->Real,
                  Symmetries->If[type===Real, Symmetric, Hermitian],
                  Bundle->name, Quiet->quiet ];
     name/: Metric[name] = metric;
     (* Tensors that are metrics have additional properties. *)
     defBunMetric[metric,name];
     If[ bname =!= Automatic || cobname =!= Automatic,
         (*then*) defineBasisNames[bname,cobname,name,quiet]];
     name/: MetricType[name] = mtype;
     (* If the user requests, make this bundle Riemannian.  If we use
        default names for the curvature tensors, create them in the
        user's context so they'll be saved with the user's stuff. *)
     If[ mtype===Riemannian,
           (*then*)
           riemanntensor = If[ rm=!=Null, rm, defaultRm];
           riccitensor   = If[ rc=!=Null, rc, defaultRc];
           scalarcurv    = If[ sc=!=Null, sc, defaultSc];
           makeRiemannian[ name, riemanntensor, riccitensor, scalarcurv, 
                           rmc, quiet ] ]; 
     (* Now protect the bundle name against accidental modification *)
     Protect[name];
     DeclareBundle[name, TorsionFree -> tfree,     
                         FlatConnection -> flat,
                         ParallelFrame -> pframe,  
                         OrthonormalFrame -> oframe,
                         CommutingFrame -> cframe, 
                         PositiveDefinite -> pdef];
     (* Tell the user what we've done. *)
     If[!quiet,
       Print["Bundle ",name," defined."];
       Print["   Metric = ",Metric[name],
                 "   Dimension = ",Dimension[name],
                 "   Indices = ",BundleIndices[name]];
       Print["   Bundle Type = ",Type[name],
         "   Metric Type = ",mtype];
       Print["   Tangent Bundle = ", TangentBundle[name]];
       Print["   Connection is ",If[flat, "flat, ", ""],
                If[!tfree, "not ", ""],"torsion free."];
       ];
     Return[name]
  ]];

(* The following function defines the extra rules needed by metric 
   tensors. *)

defBunMetric[Tensor[metric_,{},{}],bundlename_] := (
  Unprotect[metric];
  metric/: TensorMetricQ[metric] = True;
  (** Define the rule
              j                j
      metric    := Kronecker
            i               i      **)
  metric/: HoldPattern[Tensor[metric,{L[i_],U[j_]},{}]]
           := Tensor[Kronecker,{L[i],U[j]},{}];
  metric/: HoldPattern[Tensor[metric,{U[i_],L[j_]},{}]]
           := Tensor[Kronecker,{U[i],L[j]},{}];
  (** Always assume the metric is parallel. **)
  metric/: HoldPattern[Tensor[metric,{i_,j_},{__}]] := 0;
  (** If the frame is orthonormal, derivatives of the metric vanish **)
  If[ OrthonormalFrame[bundlename],
      metric/: HoldPattern[Tensor[metric,{_,_},{__}]] := 0;
      (** For 1-dimensional bundles with orthonormal frames, 
          the metric coefficients are all 1 **)
      If[ Dimension[bundlename] === 1,
          metric/: HoldPattern[Tensor[metric,{_,_},{}]] := 1;
        ]
    ];
  Protect[metric];
  );

defineBasisNames[bname_,cobname_,bundle_,quiet_] := (
  Unprotect[{Basis,bname}];
  bname/: HoldPattern[bname[L[i_]]] := 
            Tensor[Basis,{L[i]},{}];
  Basis/: HoldPattern[TensorFormat[Basis,{L[i_]},{}]] := 
            TensorFormat[bname,{L[i]},{}] /;
            MemberQ[{bundle,Conjugate[bundle]},Bundle[i]];
  Basis/: HoldPattern[TensorTeXFormat[Basis,{L[i_]},{}]] := 
            TensorTeXFormat[bname,{L[i]},{}] /;
            MemberQ[{bundle,Conjugate[bundle]},Bundle[i]];
  If[ cobname === Automatic,
      (*then*) bname/: HoldPattern[bname[U[i_]]] := 
                         Tensor[Basis,{U[i]},{}];
               Basis/: HoldPattern[TensorFormat[Basis,{U[i_]},{}]] := 
                         TensorFormat[bname,{U[i]},{}] /;
                         MemberQ[{bundle,Conjugate[bundle]},Bundle[i]];
               Basis/: HoldPattern[TensorTeXFormat[Basis,{U[i_]},{}]] := 
                         TensorTeXFormat[bname,{U[i]},{}] /;
                         MemberQ[{bundle,Conjugate[bundle]},Bundle[i]],
      (*else*) Unprotect[cobname];
               cobname/: HoldPattern[cobname[U[i_]]] :=
                         Tensor[Basis,{U[i]},{}];
               Basis/: HoldPattern[TensorFormat[Basis,{U[i_]},{}]] := 
                         TensorFormat[cobname,{U[i]},{}] /;
                         MemberQ[{bundle,Conjugate[bundle]},Bundle[i]];
               Basis/: HoldPattern[TensorTeXFormat[Basis,{U[i_]},{}]] := 
                         TensorTeXFormat[cobname,{U[i]},{}] /;
                         MemberQ[{bundle,Conjugate[bundle]},Bundle[i]];
               Protect[cobname];
    ];
  Protect[{Basis,bname}];
  );

(******************* DeclareBundle ****************************)

Options[DeclareBundle] := {FlatConnection   -> Null,
                           ParallelFrame    -> Null,
                           OrthonormalFrame -> Null,
                           PositiveDefinite -> Null,
                           CommutingFrame   -> Null,
                           TorsionFree      -> Null};

DeclareBundle[ b_, opts___ ] := 
  Module[ {flat,tfree,pframe,oframe,cframe,pdef},
  {flat,tfree,pframe,oframe,cframe,pdef} =
    {FlatConnection,TorsionFree,ParallelFrame,
     OrthonormalFrame,CommutingFrame,PositiveDefinite} /.
    {opts} /. Options[DeclareBundle];
  Which[
   flat =!= Null && !tfQ[flat],
     Message[DefineBundle::invalid,flat,"value for FlatConnection"],
   tfree =!= Null && !tfQ[tfree],
     Message[DefineBundle::invalid,tfree,"value for TorsionFree"],
   pframe =!= Null && !tfQ[pframe],
     Message[DefineBundle::invalid,pframe,"value for ParallelFrame"],
   oframe =!= Null && !tfQ[oframe],
     Message[DefineBundle::invalid,oframe,"value for OrthonormalFrame"],
   cframe =!= Null && !tfQ[cframe],
     Message[DefineBundle::invalid,cframe,"value for CommutingFrame"],
   pdef =!= Null && !tfQ[pdef],
     Message[DefineBundle::invalid,pdef,"value for PositiveDefinite"],
   (* Now all the arguments are OK.  Let's change the options. *)
   True,
     Unprotect[b];
     If[ tfree  =!= Null,  b/: TorsionFree[b] = tfree ];
     If[ pframe ||
         flat   =!= Null,  b/: FlatConnection[b] = flat ];
     If[ pframe =!= Null,  b/: ParallelFrame[b] = pframe ];
     If[ oframe =!= Null,  b/: OrthonormalFrame[b] = oframe ];
     If[ cframe =!= Null,  b/: CommutingFrame[b] = cframe ];
     If[ pdef   =!= Null,  b/: PositiveDefinite[b] = pdef ];
     Protect[b];
     Return[b]
  ] ];
  
(******************* UndefineBundle ****************************)

Attributes[UndefineBundle] = {Listable};

Options[UndefineBundle] := {Quiet -> $Quiet};

UndefineBundle[bundle_,opts___] := Module[{quiet},
  If[!optionsOK[UndefineBundle,opts],Return[]];
  quiet = Quiet /. {opts} /. Options[UndefineBundle];
  If[!BundleQ[bundle],
     Message[UndefineBundle::invalid,bundle,"bundle name"];
     Return[]
    ];
  UndefineIndex[BundleIndices[bundle],Quiet->quiet];
  UndefineTensor[Metric[bundle],Quiet->quiet];
  If[ MetricType[bundle] === Riemannian,
      UndefineTensor[RiemannTensor[bundle]];
      UndefineTensor[RicciTensor[bundle]];
      UndefineTensor[ScalarCurv[bundle]] ];
  If[ !quiet , Print["Undefining bundle ",bundle]];
  Unprotect[bundle];
  ClearAll[bundle];
  ];

(*************************** Constant *****************************)

(********************* ConstantQ ***********************************)

(* This function simply checks whether there are any tensors
   in the expression.  The reason for not using the Constant
   attribute to distinguish constants is that a symbol that
   is a constant in the "space" variables (i.e. with respect to
   covariant differentiation) may be a parameter which should
   act as a variable when applying D or Dt. *)

ConstantQ[x_] := FreeQ[x,Tensor[_,_,_]];

(******************** DefineConstant *******************************)

SetAttributes[DefineConstant,{Listable}];

Options[DefineConstant] := {Type -> Real,
                            Quiet -> $Quiet};

contypeQ[t_List] := 
  (t==={Complex}) ||
  (t==={Imaginary}) ||
  (ContainsQ[ {Real,Positive,Negative,NonPositive,
               NonNegative,Integer,Even,Odd}, t ] &&
   !ContainsQ[ t, {Positive,Negative} ] &&
   !ContainsQ[ t, {Positive,NonPositive} ] &&
   !ContainsQ[ t, {Negative,NonNegative} ] &&
   !ContainsQ[ t, {Odd,Even} ]);
contypeQ[t_] := contypeQ[{t}];

DefineConstant[c_, opts___] :=
  Module[{type,quiet},
    If[!optionsOK[DefineConstant,opts],Return[]];
    {type,quiet} = {Type,Quiet} /. {opts} /. Options[DefineConstant];
    Which[
      !newsymQ[c],   
        Message[DefineConstant::invalid,HoldForm[c],"constant name"],
      !contypeQ[type],
        Message[DefineConstant::invalid,type,"constant type"],
      True,
        DeclareConstant[c,Type->type];
        c/: c[i:(_L|_U)...] := InsertIndices[c,{i}];
        If[ !quiet, 
            Print["Constant ",c," defined.  Conjugate[",c,"] = ",
                   Conjugate[c],"."]
          ];
    ];
    Return[c];
  ];

(********************* UndefineConstant *******************************)

SetAttributes[UndefineConstant, {Listable,HoldFirst}];

Options[UndefineConstant] := {Quiet -> $Quiet};

UndefineConstant[c_,opts___] := Module[{quiet},
  If[!optionsOK[UndefineConstant,opts],Return[]];
  quiet = Quiet /. {opts} /. Options[UndefineBundle];
  Which[
    !optionsOK[UndefineConstant,opts],  
               Return[],
    !(Head[c]===Symbol && ConstantQ[c]),
               Message[UndefineConstant::invalid,HoldForm[c],"constant"],
    True,      Clear[c];
               If[!quiet,Print["Constant ",c," undefined."]]
    ]
  ];

(********************* Declare, DeclareConstant *************************)

Attributes[Declare] = {Listable};

Options[Declare] := Union[ Options[DeclareTensor],
                           Options[DeclareBundle],
                           Options[DeclareMathFunction],
                           Options[DeclareIndex],
                           Options[DeclareConstant] ];

Declare[x_, opts___] := Module[{sub},
  Switch[x,
    Tensor[t_,{},{}],       sub = DeclareTensor,
    b_Symbol?BundleQ,       sub = DeclareBundle,
    i_Symbol?IndexQ,        sub = DeclareIndex,
    f_Symbol?mathFunctionQ, sub = DeclareMathFunction,
    c_Symbol?ConstantQ,     sub = DeclareConstant,
    _,                      Message[Declare::invalid,x,
                              "tensor, bundle, index, math function, or constant"];
                            Return[];
    ];
  If[ !optionsOK[sub,opts], Return[]];
  sub[x,opts]
  ];

Options[DeclareConstant] := {Type -> Null};

DeclareConstant[ c_, opts___ ] := Module[{type},
  If[MemberQ[optionNames[{opts}], Type],
    type = Flatten[{Type /. {opts}}];
    If[ !contypeQ[type], 
        (*then*) Message[Declare::invalid, type, "constant type"];
                 Return[]
      ];
    If[ ValueQ[Sign[c]],        c/: Sign[c]        =.];
    If[ ValueQ[Positive[c]],    c/: Positive[c]    =.];
    If[ ValueQ[NonNegative[c]], c/: NonNegative[c] =.];
    If[ ValueQ[NonPositive[c]], c/: NonPositive[c] =.];
    If[ ValueQ[Negative[c]],    c/: Negative[c]    =.];
    If[ ValueQ[Conjugate[c]],   c/: Conjugate[c]   =.];
    If[ EvenQ[c],               c/: EvenQ[c]       =.];
    If[ OddQ[c],                c/: OddQ[c]        =.];
    If[ IntegerQ[c],            c/: IntegerQ[c]    =.];
    If[ MemberQ[type, Real|Positive|Negative|
                      NonNegative|NonPositive|
                      Integer|Even|Odd],       c/: Conjugate[c]   = c ];
    If[ MemberQ[type, Imaginary],              c/: Conjugate[c]   = -c ];
    If[ MemberQ[type, Positive],               c/: Sign[c]        = 1;
                                               c/: Positive[c]    = True;
                                               c/: NonNegative[c] = True ];
    If[ MemberQ[type, NonNegative],            c/: NonNegative[c] = True;
                                               c/: Negative[c]    = False ];
    If[ MemberQ[type, NonPositive],            c/: Positive[c]    = False;
                                               c/: NonPositive[c] = True ];
    If[ MemberQ[type, Negative],               c/: Sign[c]        = -1;
                                               c/: Positive[c]    = False;
                                               c/: NonNegative[c] = False ];
    If[ MemberQ[type, Integer|Even|Odd],       c/: IntegerQ[c]    = True ];
    If[ MemberQ[type, Even],                   c/: EvenQ[c] = True ];
    If[ MemberQ[type, Odd],                    c/: OddQ[c]  = True ];
  ];
  Return[c] ];


(*************************** DefineRelation *****************************)

(************ DefineRelation, UndefineRelation, & DefineRule ****************)

(* These commands are both implemented by the internal subroutine defineR.
   The exact function to perform is indicated by the argument to defineR:

   rhs = Null:       UndefineRelation
   rulename = Null:  DefineRelation
   rulename != Null: DefineRule
*)

SetAttributes[{DefineRelation,UndefineRelation,DefineRule},HoldAll];

Options[DefineRelation] := {Quiet -> $Quiet};

Options[UndefineRelation] := {Quiet -> $Quiet};

Options[DefineRule] := {Quiet   -> $Quiet,
                        NewRule -> False};

DefineRelation[]    := Message[DefineRelation::argm,"DefineRelation",0,2];
DefineRelation[_]   := Message[DefineRelation::argmu,"DefineRelation",2];
DefineRule[]        := Message[DefineRule::argm, "DefineRule",0,3];
DefineRule[_]       := Message[DefineRule::argmu,"DefineRule",3];
DefineRule[_,_]     := Message[DefineRule::argm,"DefineRule",2,3];
UndefineRelation[]  := Message[UndefineRelation::argm,"UndefineRelation",0,1];

UndefineRelation[t_,opts___] := 
  defineR[Hold[t],Hold[Null],Hold[Null],Hold[Null],UndefineRelation,opts];

DefineRelation[lhs_,rhs_,opts___Rule] := 
  defineR[Hold[lhs],Hold[rhs],Hold[Null],Hold[Null],DefineRelation,opts];

DefineRelation[lhs_, rhs_, cond_, opts___] := 
  defineR[Hold[lhs],Hold[rhs],Hold[cond],Hold[Null],DefineRelation,opts];

DefineRule[rulename_,lhs_,rhs_,opts___Rule] := 
  defineR[Hold[lhs],Hold[rhs],Hold[Null],Hold[rulename],DefineRule,opts];

DefineRule[rulename_, lhs_, rhs_, cond_, opts___] := 
  defineR[Hold[lhs],Hold[rhs],Hold[cond],Hold[rulename],DefineRule,opts];

defineR[Hold[lhs_],Hold[rhs_],Hold[cond_],Hold[rulename_],
        fcnname_,opts___] := Module[{newrule, quiet, holdlhs,
                                    singletensorQ, lhsdums, lhsfrees,
                                    lhsones, lhsonesloc, newlhsones,
                                    newlhsonenames, newlhsonepatterns,
                                    lhsdumnames, pairlhsdumnames,
                                    lhsdumrules, lhsfreerules,
                                    rhsfreerules, lhsindexlist,
                                    bundles, lhsnames, bunpairs,
                                    holdcond, complexnames},

  If[!optionsOK[fcnname,opts], Return[]];

  {newrule,quiet} = {NewRule,Quiet} /. {opts} /. Options[fcnname];

  (* Check that rulename is valid--either a new symbol or a list of rules *)

  If[!MatchQ[Hold[rulename],Hold[_Symbol]] ||
     (ValueQ[rulename] && !MatchQ[rulename,_List]),
     Message[Ricci::invalid, rulename, "rule name"];
     Return[]];

  (* Convert input tensors in lhs to internal form; if any are in
     internal form to begin with, first put them back in input form to
     avoid infinite loops. *)

  holdlhs = Hold[lhs] //. HoldPattern[Tensor[t_,{i___},{j___}]] :> t[i][j]
                      /. {t_Symbol?TensorQ :> Tensor[t,{},{}]}
                      //. HoldPattern[Inverse[x_][i__]] :> 
                            Tensor[Inverse[x],{i},{}]
                      //. HoldPattern[Conjugate[Tensor[Inverse[x_],{},{}]]] :>
                            Evaluate[Tensor[Inverse[Conjugate[x]],{},{}]]
                      //. HoldPattern[Conjugate[Tensor[t_[Bar],{},{}]]] :>
                            Tensor[t,{},{}]
                      //. HoldPattern[Conjugate[Tensor[t_,i_,j_]]] :>
                            Tensor[t[Bar],{},{}]
                      //.{HoldPattern[Tensor[t_,{i___},{j__}][(k___)?IndexQ]] :> 
                            Tensor[t,{i},{j,k}],
                          (Tensor[t_,{i___},{}])?slotsfullQ [(j___)?IndexQ] :> 
                            Tensor[t,{i},{j}],
                          HoldPattern[Tensor[t_,{i___},{}][(j___)?IndexQ]] :> 
                            Tensor[t,{i,j},{}]}
                      /. barrule;

  (* Require that lhs is a single tensor (with or without indices)
     for DefineRelation & UndefineRelation *)

  singletensorQ = MatchQ[holdlhs,(Hold[Tensor[_,_,_]])];

  If[ (fcnname === DefineRelation || fcnname === UndefineRelation) &&
    !singletensorQ,
       (*then*) Print["ERROR: first argument to ",fcnname,
                      " must be a single tensor"];
                Return[]];

  (* OK, arguments look fine.  Let's proceed. *)

  (* NewRule -> True means erase previously existing rules by the same name *)

  If[newrule, rulename = .];

  (* BUILD LEFT-HAND SIDE.  Begin transformations of indices in lhs.
     We separate lhs indices into three groups: dummy indices, free
     indices, and one-dimensional indices.  First get separate lists
     of the three categories. *)

  lhsdums = GetDummyIndices[ holdlhs ];
  lhsfrees = GetFreeIndices[ holdlhs ];
  lhsones = Select[ lhsfrees, Dimension[Bundle[#]] === 1 & ];
  lhsfrees = Complement[ lhsfrees, lhsones ];

  (* Deal with the one-dimensional indices first.  Find exactly where
     they occur. *)

  lhsonesloc = Position[ holdlhs, Alternatives @@ lhsones ];
  lhsones = PartAt[holdlhs,#]& /@ lhsonesloc;

  (* Give new names to the one-dimensional indices, since when the
     relation is invoked, different occurrences of the index may
     appear at different altitudes. *)

  newlhsones = newdummyname /@ lhsones;
  newlhsonenames = First /@ newlhsones /. x_[Bar] :> x; 
  newlhsonepatterns = Pattern[#,Blank[]]& /@ newlhsonenames;

  (* Replace all the one-dimensional indices on lhs by their new names. *)

  holdlhs = Fold[ ReplacePart[ #1, #2[[1]], #2[[2]] ]&,
                  holdlhs,
                  Transpose[{newlhsonepatterns,lhsonesloc}] ];

  (* Next deal with dummy pairs on lhs.  We must create a rule that
     will recognize either "L[i],U[i]" or "U[i],L[i]".  This is done
     by using a pattern of the form "i_,pairi_" on the lhs, and then
     checking that pairi===Pair[i] before applying the rule.  The new
     list "lhsdums" will be of the form {i, j, ... } *)

  lhsdumnames = First /@ lhsdums /. x_[Bar] :> x;

  (* For lhsdumnames = {i, j, ... } created above, create a new list of
     names of the form newdums = {pairi, pairj, ... }. *)

  pairlhsdumnames = newpairname /@ lhsdumnames;

  (* dumpairs = { {i,pairi}, {j,pairj}, ... }. *)

  dumpairs = Transpose[{lhsdumnames,pairlhsdumnames}];

  (* "lhsdumrules" is a list of rules that change L[i], U[i] to i_,
     pairi_.  The condition built by buildcond will make sure that
     these indices form a matched pair before applying the rule. *)

  lhsdumrules = Flatten[ Replace[ #, {a_,b_} :> 
                            {(L[a]|L[a[Bar]]) :> Pattern[#,Blank[]]& [a],
                             (U[a]|U[a[Bar]]) :> Pattern[#,Blank[]]& [b]} ]& /@
                         dumpairs ];

  (* Now apply the rules to lhs and change dummies to patterns. *)

  holdlhs = holdlhs /. lhsdumrules;

  (* Finally deal with free indices on lhs.  "lhsfreerules" is a list
     of rules that will convert free lhs indices of the form L[i],
     U[i], L[i[Bar]], or U[i[Bar]] to i_.  This is how we get the rule
     to be applied regardless of whether the free indices are up or
     down. *)

  lhsfreerules = Cases[ lhsfrees /. { 
    x:((L[i_Symbol])|(U[i_Symbol])|(L[i_Symbol[Bar]])|(U[i_Symbol[Bar]])) :> 
      (x :> Pattern[#,Blank[]]& [i])},
    __RuleDelayed];

  (* Apply the rules. *)

  holdlhs = holdlhs /. lhsfreerules;

  (* Finally, if the left-hand side is a single tensor with or without
     indices, insert dummy variables ind1___ and ind2___ into the
     left-hand side as place holders for additional indices that might
     be inserted by the user. *)

  holdlhs = Replace[holdlhs, Hold[Tensor[t_,{i___},{j___}]] :> 
                             Hold[Tensor[t,{i,ind1___},{j,ind2___}]]];

  (* BUILD THE RHS.  First we replace all dummy indices on rhs by new
     dummy names, to prevent conflicts with index names on the LHS or
     those used when the rule is invoked; we also hold the RHS to
     prevent evaluation, and convert LB[_] and UB[_] to L[_[Bar]] and
     U[_[Bar]]. *)

  holdrhs = NewDummy[Hold[rhs] /. barrule];

  (* Now build a set of rules for transforming the rhs free indices. *)

  rhsfreerules = Cases[ lhsfrees /. { 
    x:((L[i_Symbol])|(U[i_Symbol])|(L[i_Symbol[Bar]])|(U[i_Symbol[Bar]])) :> 
      (x :> i)},
    __RuleDelayed ];

  (* Apply the rules. *)

  holdrhs = holdrhs /. rhsfreerules;

  (* If there were one-dimensional indices on the lhs, insert metric
     factors into rhs to compensate if the user puts them at different
     altitudes. *)

  holdrhs = Fold[ insertmetricfactor, 
                  holdrhs, 
                  Transpose[{lhsones,newlhsonenames}] ] ;

  (* Lastly, we enclose rhs in a "Module" so that dummy index names
     will be newly created every time the relation is invoked; this is
     to avoid conflict with dummy names such as "c2" in the user's
     call.  Begin by getting a list of dummy names to pass to
     "buildmodule". *)

  allindices = GetIndices[holdrhs];
  rhsdumnames = Union [ First /@ Intersection[
            (Pair /@ Cases[allindices,U[i_]/;Dimension[Bundle[i]]=!=1]), 
                     Cases[allindices,L[i_]/;Dimension[Bundle[i]]=!=1]
          ] /. x_[Bar] :> x ];
  holdrhs = buildmodule[holdrhs,rhsdumnames];

  (* BUILD THE CONDITION.  "bundles" is just a list of the bundle
     names associated with all the indices in lhs.  This will be Null
     if there's no bundle.  Similarly, "lhsnames" is the list of names
     of all indices that appear in lhs.  *)

  lhsindexlist = Flatten[Join[lhsfrees,lhsdums,newlhsones]];
  bundles = Bundle /@ Lower /@ lhsindexlist;
  lhsnames = (First /@ lhsindexlist) /. x_[Bar] :> x;

  (* bunpairs = { {i, buni}, {j, bunj}, ... }.  (Discard any entries
     for which the bundle name is Null.)  This will be used to
     construct the condition for applying the rule: the rule will get
     applied only if Bundle[i] = buni, etc.  (This has to be modified
     somewhat for complex indices.  See "buildcond" below.) *)

  bunpairs = Union[Cases[Transpose[{lhsnames,bundles}],
                         {_Symbol,_?(#=!=Null&)}]];

  (* Now build the condition for applying the rule. *)  

  holdcond = buildcond[ Hold[cond], bunpairs, dumpairs ];
  holdcond = holdcond /. rhsfreerules;

  (* DEFINE THE RELATION.  "defOne" defines one rule or relation, with
     appropriate side conditions. *)

  defOne[ holdlhs, holdrhs, holdcond, Hold[rulename]];

  (* See if we have to define the conjugate relation.  We do this only
     if (a) lhs is a single tensor, and (b) either the tensor name or
     one or more indices are complex (cxQ[holdlhs] == True). *)

  (* "complexnames" is a list of lhs index pattern names that are
     associated with complex bundles.  These
     are candidates for conjugation when the conjugate rule is
     generated. *)

  complexnames = Select[ Join[lhsnames, pairlhsdumnames],
                           Type[Bundle[#]]===Complex&
                       ];

  If[MatchQ[holdlhs,Hold[Tensor[_,_,_]]] && cxQ[holdlhs],
   (*then*) defOne[holdlhs /. conjTensorRule /. conjIndexRule,
                   holdrhs /. (conjrule /@ complexnames) /.
                              Hold[x_] :> Hold[Conjugate[x]] /; 
                                !(Hold[x]===Hold[Null]),
                   holdcond /. conjTensorRule /. conjIndexRule /. 
                               conjBundleRule,
                   Hold[rulename]];
    ];

  (* Tell the user what we've done. *)

  If[!quiet, Which[
    Hold[rhs]===Hold[Null],  
      Print["Relation undefined."],
    Hold[rulename]===Hold[Null],
      Print["Relation defined."],
    True,
      Print["Rule defined."]
    ]];
  ];

(************** buildcond ***********************************************)

(* This builds a condition for deciding whether to apply the rule.
   This consists of three parts: 
   (1) the user's condition, if specified; 
   (2) the index/bundle pairs, to make sure the user's indices are associated
       with the correct bundles;
   (3) the dummy pairs {i,pairi}, to make sure they form a paired set. 

   Checking the index/bundle pairs is complicated by the presence of complex
   indices.  Basically, if a is an index associated with the complex bundle E, 
   and the user specifies L[a] in lhs, we want the rule to be applied whenever 
   the supplied index is a lower E index or an upper Conjugate[E] index.  The 
   way to do this is to check that the LOWERED version of the index is associated
   with E.
*)

buildcond[Hold[Null], {}, {}] := Hold[Null];

buildcond[Hold[Null], {}, {dumpairs__}] :=
  Hold[And @@ (PairQ @@ # &) /@ {dumpairs}];

buildcond[Hold[Null], {bunpairs__}, {}] := 
  Hold[And @@ (Bundle @ Lower[ #[[1]] ] === #[[2]] &) /@ {bunpairs}];

buildcond[Hold[Null], {bunpairs__}, {dumpairs__}] :=
  Hold[And @@ (Bundle @ Lower[ #[[1]] ] === #[[2]] &) /@ {bunpairs} &&
       And @@ (PairQ @@ # &) /@ {dumpairs}];

buildcond[Hold[cond_], {}, {}] := Hold[cond] /. barrule;

buildcond[Hold[cond_], {}, {dumpairs__}] :=
  Hold[cond && 
       And @@ (PairQ @@ # &) /@ {dumpairs}] /. barrule;

buildcond[Hold[cond_], {bunpairs__}, {}] := 
  Hold[cond && 
       And @@ (Bundle @ Lower[ #[[1]] ] === #[[2]] &) /@ {bunpairs}] /. barrule;

buildcond[Hold[cond_], {bunpairs__}, {dumpairs__}] :=
  Hold[cond && 
       And @@ (Bundle @ Lower[ #[[1]] ] === #[[2]] &) /@ {bunpairs} &&
       And @@ (PairQ @@ # &) /@ {dumpairs}] /. barrule;

(*********************** defOne **************************************)

(* For DefineRelation and UndefineRelation, the case of a barred tensor
   on the lhs has to be handled separately, since we have to extract
   the tensor name as a tag.
*)

(* If rhs is Null, then undefine relation. *)

defOne[Hold[Tensor[t_[Bar], i_, j_]],Hold[Null],_,_] := (
  Unprotect[t];
  t/: Tensor[t[Bar], i, j]=.;
  Protect[t]
  );

defOne[Hold[Tensor[t_, i_, j_]],Hold[Null],_,_] := (
  Unprotect[t];
  t/: Tensor[t, i, j]=.;
  Protect[t]
  );

(* If rulename is Null, then define relation. *)

(* Case 1: no condition. *)

defOne[Hold[Tensor[t_[Bar], i_, j_]], Hold[rhs_], Hold[Null], 
    Hold[Null] ] := (
  Unprotect[t];
  t/: Tensor[t[Bar], i, j] := evalRelation[rhs,{ind1},{ind2}];
  Protect[t]
  );

defOne[Hold[Tensor[t_, i_, j_]], Hold[rhs_], Hold[Null], 
    Hold[Null] ] := (
  Unprotect[t];
  t/: Tensor[t, i, j] := evalRelation[rhs,{ind1},{ind2}];
  Protect[t]
  );

(* Case 2: condition to check. *)

defOne[Hold[Tensor[t_, i_, j_]], Hold[rhs_], Hold[cond_], Hold[Null] ] := (
  Unprotect[t];
  t/: Tensor[t, i, j] := evalRelation[rhs,{ind1},{ind2}] /; cond;
  Protect[t]
  );

defOne[Hold[Tensor[t_[Bar], i_, j_]], Hold[rhs_], Hold[cond_], Hold[Null]
    ] := (
  Unprotect[t];
  t/: Tensor[t[Bar], i, j] := evalRelation[rhs,{ind1},{ind2}] /; cond;
  Protect[t]
  );

(* If rulename is not Null, then defining a rule. *)

(* Case A: lhs is a single tensor.  Same 2 cases as for DefineRelation. *)

defOne[Hold[Tensor[t_,i_,j_]], Hold[rhs_], Hold[Null], Hold[rulename_] ] := 
  rulename = AppendRule[rulename, 
             HoldPattern[Tensor[t, i, j]] :> evalRelation[rhs,{ind1},{ind2}]];

defOne[Hold[Tensor[t_,i_,j_]], Hold[rhs_], Hold[cond_], Hold[rulename_] ] := 
  rulename = AppendRule[rulename,
    HoldPattern[Tensor[t, i, j]] :> evalRelation[rhs,{ind1},{ind2}] /; cond];

(* Case B: general left-hand side. Again two cases. *)

defOne[Hold[lhs_], Hold[rhs_], Hold[Null], Hold[rulename_] ] := 
  rulename = AppendRule[rulename, HoldPattern[lhs] :> NewDummy[rhs]];

defOne[Hold[lhs_], Hold[rhs_], Hold[cond_], Hold[rulename_] ] := 
  rulename = AppendRule[rulename, HoldPattern[lhs] :> NewDummy[rhs] /; cond];

defOne[___] := Print["ERROR: Internal error in DefineRule."];

(*************** evalRelation ******************************************)

(* Called when the right-hand side is evaluated, to decide whether 
   indices need to be inserted or not. *)

evalRelation[rhs_,{i___},{j___}] := 
  If[{i,j}==={}, 
        (*then*)  NewDummy[rhs],
        (*else*)  CovD[InsertIndices[NewDummy[rhs],{i}],{j}]
    ];

(****************** Misc. subroutines for defineR **********************)

newdummyname[ L[i_[Bar]] ] := L[ Unique[ToString[i]] [Bar] ];
newdummyname[ U[i_[Bar]] ] := U[ Unique[ToString[i]] [Bar] ];
newdummyname[ L[i_] ]      := L[ Unique[ToString[i]] ];
newdummyname[ U[i_] ]      := U[ Unique[ToString[i]] ];

insertmetricfactor[ exp_, {oldindex_,newindex_} ] := 
  Replace[ exp, Hold[y_] :> 
    Hold[ y * Metric[Bundle[oldindex]] [Pair[oldindex], newindex] ] ];

buildmodule[ Hold[exp_], {} ] := Hold[exp];

buildmodule[ Hold[exp_], dummies_ ] :=
    Hold[ Module[ dummies, exp ] ];    

conjTensorRule = {
   Tensor[t_[Bar],{i___},{j___}] :> 
     Tensor[t, {i}, {j}];
   Tensor[t_,{i___},{j___}] :> 
     Tensor[t[Bar], {i}, {j}] /; (Type[t])==={Complex}};

conjIndexRule = {
   L[i_[Bar]] :> L[i],
   U[i_[Bar]] :> U[i],
   L[i_] :> L[i[Bar]] /; Type[Bundle[i]]===Complex,
   U[i_] :> U[i[Bar]] /; Type[Bundle[i]]===Complex};

conjBundleRule = {
   Conjugate[b_Symbol?BundleQ] :> b,
   b_Symbol?BundleQ      :> Conjugate[b]};

conjrule[i_] := (i -> Conjugate[i]);

pname[HoldPattern[x_Pattern[Bar]]] := x[[1]];
pname[HoldPattern[_[x_Pattern[Bar]]]] := x[[1]];
pname[HoldPattern[x_Pattern]] := x[[1]];
pname[HoldPattern[_[x_Pattern]]] := x[[1]];
pname[_] := Null;

cxQ[Hold[Tensor[t_,{i___},{j___}]]] := 
  Type[t]==={Complex} ||
  Or @@ (Type[#]===Complex&) /@ 
    Bundle /@ pname /@ {i,j};

AppendRule[{r___}, rule_] := {r,rule};
AppendRule[r_, rule_] := {rule};

barrule = {HoldPattern[LB[x_]] :> L[x[Bar]], HoldPattern[UB[x_]] :> U[x[Bar]]};

newpairname[x_] := Module[{name,spellon,spell1on},
  spellon  = !MatchQ[ General::spell,  _$Off ];
  spell1on = !MatchQ[ General::spell1, _$Off ];
  Off[General::spell,General::spell1];
  name = ToExpression[ "pair" <> ToString[x] ];
  If[ spellon,  On[General::spell]];
  If[ spell1on, On[General::spell1]];
  Return[name]
  ];

(*************************** Deriv *****************************)

(************ Del[exp] = total covariant derivative of exp *****************)

Options[Del] := {Metric -> Automatic,
                 Conn   -> Conn};

Del[_?ConstantQ, ___Rule] := 0;
Del[f_Plus, con___Rule] := (Del[#,con]&) /@ f;

(* Leibniz rule for function multiples only *)

Del[(f_ /; Rank[f]===0) g_, con___Rule] := 
  f Del[g,con] + TensorProduct[g, Extd[f]];

HoldPattern[Del[Tensor[g_,{},{}]]] := 0 /; MetricQ[g];

HoldPattern[Del[g_Tensor,Metric->g_]] := 0;

HoldPattern[Del[___, Tensor[Basis,{i_},{}]]] := 0 /; ParallelFrame[Bundle[i]];

HoldPattern[Del[Summation[f_],con___Rule]] := Del[f,con] ;

HoldPattern[Del[(f_ /; Rank[f]===0), con___Rule]] := Extd[f];

Del[f_[g_],con___Rule] := f'[NewDummy[g]] * Del[g,con] /;
  mathFunctionQ[f];

Del[Conjugate[f_[g_]],con___Rule] := Conjugate[Del[f[g],con]] /;
  mathFunctionQ[f];

(*********** Del[v,exp] = covariant derivative in direction v *************)

HoldPattern[Del[f_ * v_, x_, con___Rule]] := f Del[v,x,con] /; Rank[f]===0;
HoldPattern[Del[v_ + w_, x_, con___Rule]] := Del[v,x,con] + Del[w,x,con];
Del[0, _, ___Rule] := 0;
Del[v_, x_Plus, con___Rule] := Del[v,#,con]& /@ x;
Del[v_, f_ * g_, con___Rule] := Del[v,f,con] g + f Del[v,g,con];

Del[v_, f_^p_, con___Rule] := 
  NewDummy[p] NewDummy[f]^(p-1) Del[v,f,con] + 
  Log[NewDummy[f]] f^NewDummy[p] Del[v,p,con];

HoldPattern[Del[v_, Wedge[x_,y__], con___Rule]] := 
  Wedge[Del[v,x,con],y] + Wedge[x,Del[v,Wedge[y],con]];

HoldPattern[Del[v_, TensorProduct[x_,y__], con___Rule]] := 
  TensorProduct[Del[v,x,con],y] + TensorProduct[x,Del[v,TensorProduct[y],con]];

HoldPattern[Del[v_, Sym[x_], con___Rule]] := Sym[Del[v,x,con]];
HoldPattern[Del[v_, Alt[x_], con___Rule]] := Alt[Del[v,x,con]];

HoldPattern[Del[_, Tensor[g_,{},{}]]] := 0 /; MetricQ[g];

HoldPattern[Del[_, Tensor[Kronecker, {_,_}, {}], ___]] := 0;

HoldPattern[Del[_, 
            Tensor[(g_ /; MetricQ[g] && 
                          OrthonormalFrame[ Bundles[g][[1,1]] ]),
                   {_,_},{}]]] := 0;

HoldPattern[Del[v_, Tensor[g_?MetricQ, {U[i_],U[j_]}, {}]]] := 
  - (g[U[i],U[#1]] g[U[j],U[#2]] Del[v, g[L[#1],L[#2]]]&) [ 
    NewBundleSymbol[Conjugate[Bundle[i]]], 
    NewBundleSymbol[Conjugate[Bundle[j]]] ];

HoldPattern[Del[v_, Tensor[Inverse[x_],{i_,j_},{}]]] :=
  - Plus @@ Flatten [ Outer @@ Join[
      {InsertIndices[NewDummy[Inverse[x]], {i,Pair[#1]}] *
       InsertIndices[NewDummy[Inverse[x]], {Pair[#2],j}] *
       Del[v, InsertIndices[x, {#1,#2}] ] & },
      NewIndex /@ Transpose[{Bundles[x],Variance[x]}] ]];

HoldPattern[Del[v_, Summation[f_],con___Rule]] := Del[v,f,con] ;

(* If basis elements commute, derivatives of functions commute too. *)

HoldPattern[ Del[ Tensor[Basis,{L[i_]},{}],
              Del[ Tensor[Basis,{L[j_]},{}], 
                   (f_ /; Rank[f]===0) ] ] ] :=
  Del[ Tensor[Basis,{L[j]},{}], Del[ Tensor[Basis,{L[i]},{}], f ]] /;
  Bundle[i]===Bundle[j] && CommutingFrame[Bundle[i]] &&
  !IndexOrderedQ[{L[i],L[j]}];

(* For the contraction operators Inner, Int, & Dot, Del[v,_] satisfies
   the product rule unless a non-default connection is specified.*)

HoldPattern[Del[v_, Inner[x_,y_]]] := 
  Inner[Del[v,x],y] + Inner[x,Del[v,y]];

HoldPattern[Del[v_, Int[x_,y_]]] := 
  Int[Del[v,x],y] + Int[x,Del[v,y]];
HoldPattern[Del[v_, Dot[x_,y__]]] :=
  Dot[ Del[v,x], y] + Dot[ x, Del[v, Dot[y] ] ];

HoldPattern[Del[v_, Inverse[x_], con___Rule ]] :=
  - Dot[ Dot[ NewDummy[Inverse[x]], 
              Del[v,x,con]], 
         NewDummy[Inverse[x]] ];

HoldPattern[Del[v_, Det[x_], con___Rule ]] :=
  Det[x] Tr[ Inverse[x] . Del[v,x,con] ] /; rank[x]===2;

Del[v_, c_?ConstantQ, con___Rule] := 0 /; Head[c] =!= Rule;

Del[v_, f_[g_], con___Rule] := f'[NewDummy[g]] * Del[v,g,con] /;
  mathFunctionQ[f];

Del[v_, Conjugate[f_[g_]], con___Rule] := Conjugate[Del[Conjugate[v], f[g], con]] /;
  mathFunctionQ[f];

Del/: Conjugate[HoldPattern[Del[x_, y___, con_Rule]]] := 
      Del[Conjugate[x],Conjugate[y],con];
Del/: Conjugate[HoldPattern[Del[x__]]] := Conjugate /@ Del[x];

(******** Insertion of indices **********)

(* Here's where the actual computation of Del is done. *)

(* First handle a special case: derivatives of Conn are left
   in the form Del[ Basis[_], Conn[_,_,_] ] to avoid "covariant
   derivatives" of Conn. *)

Del/: TCompute[ HoldPattern[ Del[ Tensor[Conn,{i_,j_,k_},{}] ] ], 
                {l_} ] :=
  Del[ Basis[l], Tensor[Conn,{i,j,k},{}] ];

(* Now the general case.  The main complication here arises when the
   tensor being differentiated has free indices, for then its total
   covariant derivative involves some connection coefficients.
   Because of this, any other function that computes covariant
   derivatives of arbitrary tensor expressions must call this one.  
   Also, this is where we compute components of covariant derivatives
   with respect to non-default connections.  *)

Del/: TCompute[ HoldPattern[Del[x_, con___Rule]], {j___,i_}] := 
  Module[{n,diff,freeindices,connection},
  freeindices = GetFreeIndices[x];
  Switch[ {con},
    {},                connection = Conn,
    {Connection -> _}, connection = con[[2]],
    {Metric -> _},     connection = LeviCivitaConnection[ con[[2]] ],
    _,                 Message[Del::invalid, con, "Del option"];
                       Return[ERROR[Del[x,con][j,i]]] ];
  diff = connection - Conn;
  CovD[ InsertIndices[x,{j}], {i} ] +
  (* Add non-default connection terms associated to previous indices *)
  If[ diff===0,
      0,
      Sum[
        Plus @@ Flatten [ Outer @@ Join[
        {InsertIndices[ x, Join[Take[{j},n-1],{#},Drop[{j},n]] ] *
           If[ Variance[x][[n]]===Covariant,
               (*then covariant index*) 
                 - InsertIndices[ diff, { {j}[[n]], Pair[#], i } ],
               (*else contravariant*)
                   InsertIndices[ diff, { Pair[#], {j}[[n]], i } ]
             ] & },
        {NewIndex [ {Bundles[x][[n]], Variance[{j}[[n]]]} ]} ] ],
      {n,Length[{j}]}]
    ] +
  (* Finally, add connection terms associated to free indices. *)
  Sum[ 
    InsertIndices[ x /. freeindices[[n]] -> #, {j} ] *
    If[IndexAltitude[#]===Upper,
      (*then contravariant index*)
        - InsertIndices[Conn, { Pair[#], freeindices[[n]], i } ],
      (*else covariant*)
        + InsertIndices[Conn, { freeindices[[n]], Pair[#], i } ]
      ] & @ 
    dumdex[ freeindices[[n]] ],
    { n, Length[freeindices] }
    ]
  ];

(* Now Del[v,x].  First do the simpler case when v is a basis vector. *)

Del/: TCompute[HoldPattern[Del[Tensor[Basis,{k_},{}], x_, con___Rule]], {j___}] := 
  InsertIndices[ Del[x,con], {j,k} ];

(* Now the general case. *)

Del/: TCompute[ HoldPattern[ Del[v_, x_, con___Rule] ], {j___}] := 
  Module[{vindices,i},
  If[!(Rank[v]===1),
    Print["ERROR: ",v," is not a 1-tensor"];
    Return[ERROR[Del[v,x,con][j]]]
    ];
  vindices = NewLowerIndex /@ Flatten[Bundles[v]];
  Return[ Sum[ InsertIndices[ v, {Pair[ vindices[[i]] ]} ] * 
               InsertIndices[ Del[x,con], {j,vindices[[i]]} ],
               {i,Length[vindices]}]
        ]
  ];

(*** Grad ******************************************************)

Options[Grad] := {Metric -> Automatic,
                 Conn   -> Conn};

Grad[f_Plus, opt___] := (Grad[#,opt]&) /@ f;
Grad[_?ConstantQ, opts___] := 0;

(* Leibniz rule for function multiples only *)

Grad[(f_ /; Rank[f]===0) g_, opt___] := 
  f Grad[g,opt] + TensorProduct[g, Grad[f,opt]];

HoldPattern[Grad[f_ ^ p_, opt___]] :=
  NewDummy[p] NewDummy[f]^(p-1) Grad[f,opt] + 
  Log[NewDummy[f]] f^NewDummy[p] Grad[p,opt] /; 
  Rank[f]===0 && Rank[p]===0;

HoldPattern[Grad[Tensor[g_,{},{}]]] := 0 /; MetricQ[g];

HoldPattern[Grad[g_Tensor,Metric->g_]] := 0;

HoldPattern[Grad[Summation[f_],opt___]] := Grad[f,opt] ;

HoldPattern[Grad[Tensor[Basis,{i_},{}]]] := 0 /; ParallelFrame[Bundle[i]];

Grad[f_[g_],opt___] := f'[NewDummy[g]] * Grad[g,opt] /;
  mathFunctionQ[f];

Grad[Conjugate[f_[g_]],opt___] := Conjugate[Grad[f[g],opt]] /;
  mathFunctionQ[f];

Grad[f_?ConstantQ] := 0;

Grad/: Conjugate[HoldPattern[Grad[x_, y___, opt_Rule]]] := 
       Grad[Conjugate[x],Conjugate[y],opt];
Grad/: Conjugate[HoldPattern[Grad[x__]]] := Conjugate /@ Grad[x];

Grad/: TCompute[ HoldPattern[Grad[x_, opt___]], {j___,i_}] := 
  Module[{met,con},
    {met, con} = {Metric, Conn} /. {opt} /. Options[Grad];
    If[ met === Automatic,
        (*then*) Return[ Del[x,opt] [j,i] ],
        (*else*) Return[ (Del[x,opt] . Inverse[met]) [j,i] ]
      ]
    ];

(*** Div ******************************************************)

HoldPattern[Div[f_Plus,opts___Rule]] := (Div[#,opts]&) /@ f;
HoldPattern[Div[(c_?ConstantQ) f_, opts___Rule]] := c Div[f,opts];

HoldPattern[Div[ (f_ /; Rank[f]===0) * x_, opts___Rule ]] :=
  f Div[x,opts] + 
    If[ Last[Variance[x]] === Contravariant,
        (*then*) x . Extd[f],
        (*else*) x . Grad[f, opts] ];

HoldPattern[Div[Summation[x_],opts___Rule]] := Div[x,opts];

HoldPattern[Div[t_,opts___Rule]] := 0 /; Rank[t]===0;

Options[Div] := {Metric -> Automatic,
                 Connection -> Conn};

Div/: TCompute[ HoldPattern[Div[t_,opts___Rule]], {i___}] := 
  Module[{bun,met,con},
  bun = UnderlyingTangentBundle[t];
  Which[
    MemberQ[bun,Null],
      Message[Bundle::error,Div[t,opts]];
      Return[ERROR[Div[t,opts][i]]],
    !optionsOK[Div, opts],
      Return[ERROR[Div[t,opts][i]]],
    {opts} === {},
      Return[Plus @@ ( Del[t,opts] [i, #, Pair[#] ] &) /@ NewLowerIndex /@ bun],
    True,
      {met,con} = {Metric,Connection} /. {opts} /. Options[Div];
      If[ con === Conn,
          con = LeviCivitaConnection[met]];
      If[ Last[Variance[t]] === Contravariant,
          (*then*) Return[ Plus @@ 
                             ( (Del[t,Connection->con] [i, Pair[#], # ] &)
                               /@ NewLowerIndex /@ bun ) ],
          (*else*) Return[ Plus @@ ( InsertIndices[
                         If[ met===Automatic,
                             Metric[Bundle[#[[1]]]],
                             Inverse[ met ]],
                         {Pair[ #[[1]] ],Pair[ #[[2]] ]}] *
                       Del[t,Connection->con] [i, #[[1]], #[[2]] ] &)
                         /@ (Table[NewLowerIndex[#],{2}]&) /@ bun ]
          ]    
    ];
  ];

Div/:  Conjugate[HoldPattern[Div[x_,opts___Rule]]] := 
           Div[Conjugate[x],opts];

(****************** Laplacian *****************************************)

Laplacian[x_,opts___] := 
  Switch[ $LaplacianConvention,
          DivGrad,            Div[ Grad[ x, opts ], opts],
          PositiveSpectrum, - Div[ Grad[ x, opts ], opts],
          _,                Message[ $LaplacianConvention::invalid, 
                                     $LaplacianConvention, 
                                     "$LaplacianConvention" ];
                            ERROR[ Laplacian[x,opts] ]
        ];

(****************** LaplaceBeltrami **************************************)

LaplaceBeltrami[x_,opts___] := 
  Extd[ ExtdStar[ x, opts ] ] + ExtdStar[ Extd[x], opts];

(************************ Extd ************************************)

MakeLinear[Extd];

HoldPattern[Extd[_Extd]] = 0;

HoldPattern[Extd[f_ * g_]] :=
  Wedge[Extd[f], g] + f Extd[g] /; Rank[f] === 0;

HoldPattern[Extd[Wedge[f_,g__]]] := 
  Wedge[Extd[f],g] + (-1)^Rank[f] Wedge[f,Extd[Wedge[g]]];

HoldPattern[Extd[f_ ^ p_]] :=
  NewDummy[p] NewDummy[f]^(p-1) Extd[f] + 
  Log[NewDummy[f]] f^NewDummy[p] Extd[p] /; 
  Rank[f]===0 && Rank[p]===0;

HoldPattern[Extd[Summation[f_]]] := Extd[f] ;

HoldPattern[Extd[Tensor[Kronecker,{_,_},{}]]] := 0;

HoldPattern[Extd[Tensor[(g_ /; MetricQ[g] && 
                           OrthonormalFrame[ Bundles[g][[1,1]] ]),
                    {_,_},{}]]] := 0;

HoldPattern[Extd[Tensor[Inverse[x_],{i_,j_},{}]]] :=
  - Plus @@ Flatten [ Outer @@ Join[
      {InsertIndices[NewDummy[Inverse[x]], {i,Pair[#1]}] *
       InsertIndices[NewDummy[Inverse[x]], {Pair[#2],j}] *
       Extd[ InsertIndices[x, {#1,#2}] ] & },
      NewIndex /@ Transpose[{Bundles[x],Variance[x]}] ]];

HoldPattern[Extd[Tensor[g_?MetricQ,{U[i_],U[j_]},{}]]] :=
  - (g[U[i],U[#1]] g[U[j],U[#2]] Extd[g[L[#1],L[#2]]]&) [ 
    NewBundleSymbol[Conjugate[Bundle[i]]], 
    NewBundleSymbol[Conjugate[Bundle[j]]] ];

HoldPattern[Extd[f_?ConstantQ]] := 0;

HoldPattern[Extd[f_[g_]]] := f'[NewDummy[g]] * Extd[g] /;
  mathFunctionQ[f];

HoldPattern[Extd[Conjugate[f_[g_]]]] := Conjugate[Extd[f[g]]] /;
  mathFunctionQ[f];

(** Here's where we insert indices. **)

Extd/: TCompute[ HoldPattern[Extd[f_]], {i__} ] := 
  Module[{rk,tanbun,dummyindices,m,n,k},
  rk = rank[f];
  If[!FormQ[f],
    Print["ERROR: ",f," is not a differential form"];
    Return[ERROR[Extd[f][i]]]
    ];
  If[rk===0,
    (*then*) 
      (* Here Del[f] must be kept Unevaluated, or it will be
         turned back into Extd[f] by the rules for Del. *)
      TCompute[ Unevaluated[Del[f]], {i} ],
    (*else*) 
      tanbun = UnderlyingTangentBundle[f];
      (* Here's the actual computation:
         First the numerical factor (depending on WedgeConvention *)
      Signature[{i}] * WedgeFactor[{1,rk}] / (rk+1) * (
      (* Next the antisymmetrized derivatives *)
      Sum[ 
        ((-1)^rk * Signature[#] * InsertIndices[Del[f],#] &) 
           @ RotateLeft[{i},n],
        {n,rk+1}] +
      (* Finally the torsion terms, if there are any free indices *)
      If[And @@ TorsionFree /@ tanbun,
        (*then*) 
          0,
        (*else*)
          dummyindices = NewLowerIndex /@ tanbun;
          Sum[
            ((1/2) * Signature[#] *
              InsertIndices[f, Join[ {dummyindices[[k]]}, Take[#,rk-1] ] ] * 
              InsertIndices[Tor, Join[ {Pair[dummyindices[[k]]]}, 
                                       Drop[#,rk-1] ] ]&
            )
              @ (Join[ RotateLeft[ Drop[#,-1], m ], Take[#,-1] ] & @
                RotateLeft[{i},n]),
            {k,Length[dummyindices]},
            {m,rk},
            {n,rk+1}
            ]
        ]
      ) 
    ]
  ];

Extd/: Conjugate[HoldPattern[Extd[x_]]] := Extd[Conjugate[x]];

(****************************** ExtdStar ***********************************)

HoldPattern[ExtdStar[f_Plus,opts___Rule]] := (ExtdStar[#,opts]&) /@ f;
HoldPattern[ExtdStar[(c_?ConstantQ) f_, opts___Rule]] := c ExtdStar[f,opts];

HoldPattern[ExtdStar[Summation[x_],opts___Rule]] := ExtdStar[x,opts];

HoldPattern[ExtdStar[t_,opts___Rule]] := 0 /; rank[t]===0;

HoldPattern[ExtdStar[f_*t_]] := f ExtdStar[t] - Int[Extd[f], t] /; rank[f]===0;
HoldPattern[ExtdStar[ExtdStar[f_,opts___Rule],opts___Rule]] := 0;

Options[ExtdStar] := {Metric -> Automatic};

ExtdStar/: TCompute[ HoldPattern[ExtdStar[f_,opt___Rule]], {i___} ] := 
  Module[{rk,tanbun,dummyindices,m,n,k,met,con},
  If[ !FormQ[f],
    Print["ERROR: ",f," is not a differential form"];
    Return[ERROR[ExtdStar[f][i]]]];
  rk = Rank[f];
  tanbun = UnderlyingTangentBundle[f];
  If[ MemberQ[bun,Null],
      Message[Bundle::error,ExtdStar[t,opt]];
      Return[ERROR[ExtdStar[t,opt][i]]]];
  If[ !optionsOK[ExtdStar, opt],
      Return[ERROR[ExtdStar[t,opt][i]]]];
  (* Here's the computation:
     First the numerical factor, depending on WedgeConvention *)
  (-1)^(rk-1) * WedgeFactor[{rk-1,1}] * 
   IntFactor[rk,rk] / IntFactor[rk-1,rk-1] *
   (* Next the derivative terms *)
   (InsertIndices[ - Div[f,opt], {i}] +
   (* Finally the torsion terms, if there are free indices *)
    If[rk === 1 || (And @@ TorsionFree /@ tanbun),
      (*then*)
        0,
      (*else*)
        dummyindices = Table[NewLowerIndex /@ tanbun, {2}];
        - (1/2) * 
        Sum[ 
         (Signature[#] *
          InsertIndices[f, 
            Join[ {dummyindices[[1,m]],dummyindices[[2,n]]},
                  Take[#,rk-2] ] ] *
          InsertIndices[Tor,
            Join[ Drop[#,rk-2], 
                  Pair /@ {dummyindices[[1,m]],dummyindices[[2,m]]} ]] &)
          @ RotateLeft[{i},k],
          {m,Length[dummyindices[[1]]]},
          {n,Length[dummyindices[[2]]]},
          {k,rk-1}]
      ]
    )
  ];

ExtdStar/:  Conjugate[HoldPattern[ExtdStar[x_]]] := ExtdStar[Conjugate[x]];

(****************** Generic structure equations ******************)

CompatibilityRule := {
  HoldPattern[Extd[Tensor[g_?MetricQ,{L[i_],L[j_]},{}]]] :>
    Tensor[Conn,{L[i],L[j]},{}] + Tensor[Conn,{L[j],L[i]},{}],
  HoldPattern[Del[v_,Tensor[g_?MetricQ,{L[i_],L[j_]},{}]]] :>
    Tensor[Conn,{L[i],L[j]},{}] . v +
    Tensor[Conn,{L[j],L[i]},{}] . v,
  HoldPattern[Extd[Tensor[g_?MetricQ,{U[i_],U[j_]},{}]]] :> 
    - Tensor[Conn,{U[i],U[j]},{}] - Tensor[Conn,{U[j],U[i]},{}],
  HoldPattern[Del[v_,Tensor[g_?MetricQ,{U[i_],U[j_]},{}]]] :>
    - Tensor[Conn,{U[i],U[j]},{}] . v -
      Tensor[Conn,{U[j],U[i]},{}] . v  
  };

FirstStructureRule := {HoldPattern[Extd[Tensor[Basis,{U[i_]},{}]]] :>
  (WedgeFactor[2]/2) Tensor[Tor,{U[i]},{}] +
  (Tensor[Basis,{Pair[#]},{}] ~Wedge~ Tensor[Conn,{#,U[i]},{}] &)
   @ NewLowerIndex[Bundle[i]] };

SecondStructureRule := {
  HoldPattern[ Extd[ Tensor[Conn,{L[i_],U[j_]},{}] ] ] :>
    (WedgeFactor[2]/2) Tensor[Curv,{L[i],U[j]},{}] +
    (Tensor[Conn,{L[i],Pair[#]},{}] ~Wedge~ Tensor[Conn,{#,U[j]},{}] &)
      @ NewLowerIndex[Bundle[i]],
  HoldPattern[ Del[ Tensor[Basis,{L[k_]},{}],
                Tensor[Conn,{L[i_],U[j_],L[l_]},{}] ] *
           Wedge[ a___, 
                  Tensor[Basis,{U[k_]},{}],
                  b___,
                  Tensor[Basis,{U[l_]},{}],
                  c___ ] ] :>
    ( (1/2) Curv[L[i],U[j],L[k],L[l]] +
    (Conn[L[i],U[#],L[k]] Conn[L[#],U[j],L[l]] &) @
        NewBundleSymbol[Bundle[i]] -
      Plus @@ (
        (Conn[L[i],U[j],L[#]] Conn[L[k],U[#],L[l]] +
          (1/2) Conn[L[i],U[j],L[#]] Tor[U[#],L[k],L[l]] &) /@
        NewBundleSymbol /@
        TangentBundle[Bundle[i]]
        ) ) *
    Wedge[ a,Tensor[Basis,{U[k]},{}],b,Tensor[Basis,{U[l]},{}],c],
  HoldPattern[ Del[ Tensor[Basis,{L[k_]},{}],
                Tensor[Conn,{L[i_],U[j_],L[l_]},{}] ] *
           Wedge[ a___, 
                  Tensor[Basis,{U[l_]},{}],
                  b___,
                  Tensor[Basis,{U[k_]},{}],
                  c___ ] ] :>
    ( (1/2) Curv[L[i],U[j],L[k],L[l]] +
    (Conn[L[i],U[#],L[k]] Conn[L[#],U[j],L[l]] &) @
        NewBundleSymbol[Bundle[i]] -
      Plus @@ (
        (Conn[L[i],U[j],L[#]] Conn[L[k],U[#],L[l]] +
          (1/2) Conn[L[i],U[j],L[#]] Tor[U[#],L[k],L[l]] &) /@
        NewBundleSymbol /@
        TangentBundle[Bundle[i]]
        ) ) *
    Wedge[ a,Tensor[Basis,{U[l]},{}],b,Tensor[Basis,{U[k]},{}],c],
  HoldPattern[ a_. * Del[ Tensor[Basis,{L[k_]},{}],
                      Tensor[Conn,{L[i_],U[j_],L[l_]},{}] ] +
           b_. * Del[ Tensor[Basis,{L[l_]},{}],
                      Tensor[Conn,{L[i_],U[j_],L[k_]},{}] ] ] :>
    b * (
      Plus @@ (
      (Conn[L[i], U[j], L[#]]*Tor[U[#], L[k], L[l]] +
       Conn[L[i], U[j], L[#]]*Conn[L[k], U[#], L[l]] - 
       Conn[L[i], U[j], L[#]]*Conn[L[l], U[#], L[k]] & ) /@
         NewBundleSymbol /@
         TangentBundle[Bundle[i]] ) + 
      (Conn[L[i], U[#], L[l]]*Conn[L[#], U[j], L[k]] - 
       Conn[L[i], U[#], L[k]]*Conn[L[#], U[j], L[l]] & ) @
         NewBundleSymbol[Bundle[i]] - 
      Curv[L[i], U[j], L[k], L[l]]
      ) /; a+b===0 
  };

StructureRules := 
  Join[CompatibilityRule,FirstStructureRule,SecondStructureRule];

(******************** Curvature for arbitrary connections *******)

Curvature/: TCompute[ HoldPattern[Curvature[cn_]],
                      {a_,b_,i_,j_}] := Module[{diff},
  diff = cn - Conn;
  TensorExpand[
  Plus @@ Flatten[ {
   InsertIndices[ Curv, {a,b,i,j} ] +
   CovD[ InsertIndices[ diff, {a, b, j} ], {i}] - 
     CovD[ InsertIndices[ diff, {a, b, i} ], {j}],
   (InsertIndices[ diff, {a, U[#], j} ] * 
       InsertIndices[ diff, {L[#], b, i} ] -
    InsertIndices[ diff, {a, U[#], i} ] * 
       InsertIndices[ diff, {L[#], b, j} ] &) /@ 
         NewBundleSymbol /@ Bundles[diff][[1]],
   (InsertIndices[ diff, {a, b, L[#]} ] * 
       Tor[U[#], i, j] &) /@
         NewBundleSymbol /@ TangentBundle[Bundle[a]]}
  ]]
  ];

Curvature/: TCompute[ HoldPattern[Curvature[conn_]], 
                      {i_,j_}] := 
  TensorExpand[
    Plus @@ Flatten[ Outer @@ Join[ {
      InsertIndices[ Curvature[conn], {i,j,L[#1],L[#2]} ] * 
        (1/WedgeFactor[2]) *
        Wedge[ Basis[U[#1]], Basis[U[#2]] ]
      &},
      Table[ NewBundleSymbol /@ UnderlyingTangentBundle[conn[L[i],L[j]]], {2} ]
  ]]];

Curvature/: TCompute[ HoldPattern[Curvature[h_]], {i___} ] := (
  Message[ Index::error, "Curvature" ];
  Return[ERROR[ InsertIndices[ Curvature[h], {i} ]]];
  ) /; Length[{i}] =!= 2 && Length[{i}] =!= 4;

(********************* Lie **************************************)

Lie[v_ + w_, t_] := Lie[v,t] + Lie[w,t];
Lie[(c_?ConstantQ) * v_, t_] := c Lie[v,t];
Lie[f_ * v_, w_] := f Lie[v,w] - Inner[Extd[f],w] * v /; 
  VectorFieldQ[v] && VectorFieldQ[w];
Lie[0, _] := 0;
Lie[v_, t_Plus] := Lie[v,#]& /@ t;
Lie[v_, f_ * t_] := Lie[v,f] * t + f * Lie[v,t];

Lie[v_, f_^p_] := 
  NewDummy[p] NewDummy[f]^(p-1) Lie[v,f] + 
  Log[NewDummy[f]] f^NewDummy[p] Lie[v,p];

HoldPattern[Lie[v_, Wedge[x_,y__]]] := 
  Wedge[Lie[v,x],y] + Wedge[x,Lie[v,Wedge[y]]];

HoldPattern[Lie[v_, TensorProduct[x_,y__]]] :=
  TensorProduct[Lie[v,x],y] + TensorProduct[x,Lie[v,TensorProduct[y]]];

HoldPattern[Lie[v_, Inverse[x_]]] := 
  - Dot[ Dot[ NewDummy[Inverse[x]], 
              Lie[v,x]], 
         NewDummy[Inverse[x]] ]; 

HoldPattern[Lie[v_, Sym[x_]]] := Sym[Lie[v,x]];
HoldPattern[Lie[v_, Alt[x_]]] := Alt[Lie[v,x]];

HoldPattern[Lie[v_, f_]] := Del[v,f] /; Rank[f]===0 && VectorFieldQ[v];
HoldPattern[Lie[v_, w_]] := - Lie[w, v] /; 
  VectorFieldQ[v] && VectorFieldQ[w] && !OrderedQ[{v,w}];
HoldPattern[Lie[v_, v_]] := 0;

HoldPattern[ Lie[ Tensor[Basis,{L[i_]},{}], Tensor[Basis,{L[j_]},{}] ] ] :=
  0 /; Bundle[i]===Bundle[j] && CommutingFrame[Bundle[i]];

Lie/: Conjugate[HoldPattern[Lie[x__]]] := Conjugate /@ Lie[x];

Lie/: TCompute[ HoldPattern[Lie[v_, t_]], {i___} ] := 
  Module[{vindices,tindices,variances,m,n},
  If[!VectorFieldQ[v],
    Print["ERROR: ",v," is not a vector field."];
    Return[ ERROR[ Lie[v,t][i] ] ] ];
  vindices = NewLowerIndex /@ UnderlyingTangentBundle[t];
  tindices = NewLowerIndex /@ UnderlyingTangentBundle[t];
  variances = Variance[t];
  InsertIndices[Del[v,t],{i}] +
    Sum[ 
      If[variances[[n]]===Covariant,
        (*then covariant index*)
          InsertIndices[t, Join[ Take[{i},n-1], {vindices[[m]]}, Drop[{i},n] ] ] *
          ( InsertIndices[ Unevaluated[Del[v]], 
                           {Pair[vindices[[m]]], {i}[[n]] } ] +
            Plus @@ ( InsertIndices[ v, {Pair[#]} ] * 
                      InsertIndices[ Tor, {Pair[vindices[[m]]], 
                                           #, 
                                           {i}[[n]]} ] & /@ tindices)
          ),
        (*else contravariant index*)
          InsertIndices[t, Join[ Take[{i},n-1], {Pair[vindices[[m]]]}, Drop[{i},n] ] ] *
          (- InsertIndices[ Unevaluated[Del[v]],
                            { {i}[[n]], vindices[[m]] } ] -
            Plus @@ ( InsertIndices[ v, {Pair[#]} ] * 
                      InsertIndices[ Tor, { {i}[[n]], 
                                           #, 
                                           vindices[[m]]} ] & /@ tindices)
          )
        ],
      {m, Length[vindices]},
      {n, Length[{i}]}
      ]
  ];

LieRule = 
  {HoldPattern[Lie[v_,t_]] :> 
      Int[v, Extd[t]] + Extd[Int[v,t]] /; FormQ[t]
  };

(****************** CovD *************************************************)

(* Covariant differentiation with no indices inserted is a no-op. *)

CovD[x_,{},opts___] := x;

(* First check that the arguments are OK. *)

Options[CovD] := {Metric -> Automatic,
                  Connection -> Conn};

CovD[x_, {i__}, opts___] := Which[
  !IndexQ[i],
    Print["ERROR: invalid index ",i," given to CovD"];
    Return[ERROR[CovD[x,{i},opts]]],
  !optionsOK[CovD,opts], 
    Return[ERROR[CovD[x,{i},opts]]],
  Head[rank[x]]=!=Integer,
    Print["ERROR: invalid tensor expression ",x];
    Return[ERROR[CovD[x,{i},opts]]],
  rank[x]=!=0,    
    Print["ERROR: Covariant derivative of invalid scalar expression ",x];
    Return[ERROR[CovD[x,{i},opts]]],
  True,           
    cd[x,{i},opts]
  ];

cd[x_,{i_,j__},opts___] := cd[ cd[x,{i},opts], {j}, opts];

cd[a_Plus,{i_},opts___] := cd[#,{i},opts]& /@ a;

cd[t_,{},opts___] := t;

cd[f_ * g_, {i_}, opts___] := 
  f*cd[g,{i},opts] + g*cd[f,{i},opts];

cd[f_^p_, {i_}, opts___] :=
  NewDummy[p] NewDummy[f]^(p-1) cd[f,{i},opts] + 
  Log[NewDummy[f]] f^NewDummy[p] cd[p,{i},opts];

cd[HoldPattern[Det[h_]], {i_}, opts___] :=
  InsertIndices[ Del[ Basis[i], Det[h] ], {} ];

cd[HoldPattern[Summation[f_]], {i_}] := 
  cd[f,{i}];

cd[ HoldPattern[ Tensor[Conn,{i_,j_,k_},{}] ], {l_} ] :=
  expandedCovD[ Tensor[Conn,{i,j,k},{}], {l} ];

(* The following special case has to be here because 
   Del[Basis[...],Conn[...]] is considered a component expression. *)

cd[ HoldPattern[ Del[ Tensor[Basis,{i_},{}], x_ ] ], {l_} ] :=
  expandedCovD[ Del[Tensor[Basis,{i},{}],x], {l} ];

cd[ HoldPattern[Tensor[t__]], {k_} ] := 
  TCompute[Tensor[t],{k}];

cd[f_[g_],{i_},opt___] := f'[NewDummy[g]] * cd[g,{i},opt] /;
  mathFunctionQ[f];

cd[Conjugate[f_[g_]],{i_},opt___] := Conjugate[cd[f[g],Conjugate/@{i},opt]]  /;
  mathFunctionQ[f];

cd[Tensor[t_,{a___},{b___}],{k_},Connection->con_] := 
  InsertIndices[ Del[ Nest[Del,t,Length[{b}]], Connection->con], 
                 {a,b,k} ];

cd[Tensor[t_,{a___},{b___}],{L[k_]},Metric->metric_] := 
  cd[Tensor[t,{a},{b}],{L[k]},Connection->LeviCivitaConnection[metric]];

cd[Tensor[t_,{a___},{b___}],{U[k_]},Metric->metric_] := 
  (InsertIndices[Inverse[metric], {U[k],#}] *
    cd[Tensor[t,{a},{b}],{Pair[#]},Connection->LeviCivitaConnection[metric]]& @
    dumdex[ U[k] ]
    );

cd[c_?ConstantQ,{i_},opts___] := 0;

(* Covariant derivative of anything else produces an error message. *)

cd[x_,{i__},opts___] := (
  Print["ERROR: Covariant derivative of invalid component expression "<>x];
  ERROR[ CovD[x,{i},opts] ]
  );

(*************************** subroutines **************************)

(** various subroutines used by the differentiation functions *)

(* dumdex[index] returns a new or upper index for Bundle[index]. *)

dumdex[L[i_]] :=
  NewLowerIndex[Bundle[i]];

dumdex[U[i_]] :=
  NewUpperIndex[Bundle[i]];

(* tandumdex[index] gives a list of new indices for the tangent
   decomposition. *)

tandumdex[L[i_]] := 
  NewLowerIndex /@ TangentBundle[Bundle[i]];

tandumdex[U[i_]] :=
  NewUpperIndex /@ TangentBundle[Bundle[i]];

(********************* Format *********************************)

SetAttributes[{TensorFormat,IndexFormat,IndexTeXFormat,TensorTeXFormat,
               BracketFormat,TimesFormat,InputTensorFormat,
               SymmetricProductFormat,TensorPowerFormat,
               holdlist,indexformQ,SubscriptedNameFormat
               },
              HoldAll];

(** All of the following formatting routines have to be very careful
    not to evaluate their arguments.  That's why there are so many
    HoldForm's and all these functions have the HoldAll attribute. **)

(**** Decide which format to use for tensors. ****)

Tensor/: HoldPattern[Format[Tensor[t_,{i___},{j___}]]] :=
           TensorFormat[t,{i},{j}] /; 
           indexformQ[i,j] && $TensorFormatting;

Tensor/: HoldPattern[Format[Tensor[t_,{i___},{j___}],TeXForm]] :=
           TensorTeXFormat[t,{i},{j}] /; 
           indexformQ[i,j] && $TensorTeXFormatting;

Tensor/: HoldPattern[Format[Tensor[t_,{i___},{j___}],InputForm]] :=
           InputTensorFormat[t,{i},{j}] /; 
           indexformQ[i,j] && $TensorFormatting;

(** indexformQ[i,j,...] does the same thing as IndexQ, but without
    evaluating its arguments. **)

indexformQ[] := True;
indexformQ[L[_Symbol]] := True;
indexformQ[U[_Symbol]] := True;
indexformQ[L[_Symbol[Bar]]] := True;
indexformQ[U[_Symbol[Bar]]] := True;
indexformQ[_] := False;
indexformQ[i_,j__] := indexformQ[i] && indexformQ[j];

(************************* Tensor Formatting *******************)

(* Inverse tensors with indices have to be handled specially. *)

HoldPattern[TensorFormat[Inverse[Tensor[t_,{},{}]],{i_,j_},{}]] := 
  SequenceForm[ "(", Format[Inverse[Tensor[t,{},{}]]], ")", 
                IndexFormat[{i,j},{}] ];

HoldPattern[TensorFormat[Inverse[x_],{i_,j_},{}]] := 
  SequenceForm[ "(", Format[Inverse[HoldForm[x]]], ")", 
                IndexFormat[{i,j},{}] ];

TensorFormat[t_,{i__},{}] := 
  SequenceForm[ TensorFormat[t], 
                IndexFormat[{i},{}] ];

TensorFormat[t_,{i___},{j__}] := 
  SequenceForm[ TensorFormat[t], 
                IndexFormat[{i},{j}] ];

TensorFormat[t_[Bar],{},{}] := 
  BarFormat[TensorFormat[t]];

TensorFormat[t_[Bar]] := 
  BarFormat[TensorFormat[t]];

TensorFormat[t_,{},{}] := TensorFormat[t];

TensorFormat[t_] := HoldForm[t];

(*********** TeX formatting of tensors ***************************)

HoldPattern[TensorTeXFormat[Inverse[Tensor[t_,{},{}]],{i_,j_},{}]] := 
  SequenceForm[ "\\tensor{(", Format[Inverse[Tensor[t,{},{}]],TeXForm], ")}{", 
                IndexTeXFormat[{i,j},{}], "}" ];

HoldPattern[TensorTeXFormat[Inverse[x_],{i_,j_},{}]] := 
  SequenceForm[ "\\tensor{(", Format[Inverse[x],TeXForm], ")}{", 
                IndexTeXFormat[{i,j},{}], "}" ];

TensorTeXFormat[t_,{},{}] := TensorTeXFormat[t];

TensorTeXFormat[t_,{i___},{j___}] :=
  SequenceForm[ "\\tensor{", TensorTeXFormat[t], "}{",
         IndexTeXFormat[{i},{j}], "}" ];

TensorTeXFormat[t_[Bar]] :=
  SequenceForm["\\overline{",TensorTeXFormat[t],"}"];

TensorTeXFormat[t_] := (TeXFormat /. TensorData[t] /. 
  {"" -> HoldForm[t]}) /; TensorQ[t];

TensorTeXFormat[t_] := Format[HoldForm[t],TeXForm];

(************************ Index Formatting ***************************)

IndexFormat[{L[i_[Bar]]},{}] :=
  SequenceForm[ ColumnForm[{"_",HoldForm[i]},Left,Below], " " ];

IndexFormat[{U[i_[Bar]]},{}] :=
  SequenceForm[ ColumnForm[{"_",HoldForm[i]," "},Left,Above], " " ];

IndexFormat[{L[i_]},{}] := Subscript[SequenceForm[HoldForm[i]," "]];

IndexFormat[{U[i_]},{}] := Superscript[SequenceForm[HoldForm[i]," "]];

IndexFormat[{},{}] := "";

IndexFormat[{i_,j__},{}] := 
  SequenceForm[  IndexFormat[{i},{}],IndexFormat[{j},{}] ];

IndexFormat[{i___},{j__}] := 
  SequenceForm[ IndexFormat[{i},{}],Subscript[";"],IndexFormat[{j},{}] ];

(*********************** Input format for tensors **********************)

InputTensorFormat[Inverse[Tensor[t_,{},{}]],{i_,j_},{}] := 
  Format[SequenceForm[ "Inverse[",
                       HoldForm[t],
                       "]",
                       BracketFormat[{i,j},{}] ],
         OutputForm];

InputTensorFormat[Inverse[x_],{i_,j_},{}] := 
  Format[
    SequenceForm[ Format[Inverse[x],InputForm],
                  BracketFormat[{i,j},{}] ],
    OutputForm];

InputTensorFormat[t_,{i__},{}] := 
  Format[
    SequenceForm[ InputTensorFormat[t,{},{}], 
                  BracketFormat[{i},{}] ],
    OutputForm];

InputTensorFormat[t_,{i___},{j__}] := 
  Format[
    SequenceForm[ InputTensorFormat[t,{},{}], 
                  BracketFormat[{i},{j}] ],
    OutputForm];

InputTensorFormat[t_[Bar],{},{}] := 
  Format[SequenceForm["Conjugate[",HoldForm[t],"]"],OutputForm];

InputTensorFormat[t_,{},{}] := Format[HoldForm[t],OutputForm];

(********************** Input format for indices *******************)

BracketFormat[{i___},{j__}] := 
  SequenceForm[ BracketFormat[{i},{}], " ", 
                BracketFormat[{j},{}] ];

BracketFormat[{},{}] := "[]";

BracketFormat[{i_},{}] := 
  SequenceForm[ "[", HoldForm[i], "]" ];

BracketFormat[{i_,j__},{}] := 
  SequenceForm[ "[", Infix[holdlist[i,j],", "], "]" ];

(********************* TeX formatting for indices ******************)

IndexTeXFormat[{U[i_[Bar]]},{}] := 
         SequenceForm["\\up{\\overline{",
                      Format[SubscriptedNameFormat[i],TeXForm],"}}"];

IndexTeXFormat[{L[i_[Bar]]},{}] :=
         SequenceForm["\\down{\\overline{",
                      Format[SubscriptedNameFormat[i],TeXForm],"}}"];

IndexTeXFormat[{U[i_]},{}] :=
         SequenceForm["\\up{",
                      Format[SubscriptedNameFormat[i],TeXForm],"}"];

IndexTeXFormat[{L[i_]},{}] := 
         SequenceForm["\\down{",
                      Format[SubscriptedNameFormat[i],TeXForm],"}"];

IndexTeXFormat[{i_,j__},{}] := 
  SequenceForm[  IndexTeXFormat[{i},{}],IndexTeXFormat[{j},{}] ];

IndexTeXFormat[{i___},{j__}] := IndexTeXFormat[{i,L[";"],j},{}];

(*************************** BarFormat **********************************)

(* This is used for printing barred symbols *)

BarFormat[x_] := 
  If[(* symbol is 1 character *) StringLength[ToString[HoldForm[x]]] == 1,
    (*then*) 
      ColumnForm[{"_",HoldForm[x]},Center,Above],
    (*else*)
      ColumnForm[{StringJoin @@ Table["_",{StringLength[ToString[HoldForm[x]]]}],
                  HoldForm[x]},Center,Above]
    ];

BarFormat[x_,TeXForm] := 
  SequenceForm["\\overline{", HoldForm[x], "}"];

(************************ SubscriptedNameFormat *************************)

(* This is used for TeX formatting of index names like i1, p37, etc.  The
   digits are typeset as a subscript to the index name. *)

SubscriptedNameFormat[i_Symbol] := Module[{chars,prefix,suffix,letterpos},
  chars = Characters[ToString[HoldForm[i]]];
  letterpos = Position[chars,_String?LetterQ,{1}];
  If[Length[letterpos]===Length[chars], Return[HoldForm[i]]];
  prefix = StringJoin @@ Take[chars,First[Last[letterpos]]];
  suffix = StringJoin @@ Drop[chars,First[Last[letterpos]]];
  Return[SequenceForm[Format[prefix,TeXForm],
               Subscript[Format[suffix,TeXForm]]]]
  ];

SubscriptedNameFormat[i_] := HoldForm[i];

(************************** Times *********************************)

Unprotect[Times];

(** Modify the order of products for printing. **)

Times/: HoldPattern[Format[Times[y___]]] := TimesFormat[{y}] /; 
  !FreeQ[Hold[{y}],Tensor[_,_,_]] && $TensorFormatting;

Times/: HoldPattern[Format[Times[y___],TeXForm]] := TimesFormat[{y},TeXForm] /; 
  !FreeQ[Hold[{y}],Tensor[_,_,_]] && $TensorTeXFormatting;

Protect[Times];

(* preclevel[m,n,p] chooses precedence levels based on the version of
   Mathematica in use: vers < 3.0 => m; 3.0 <= vers < 5 => n; vers >= 5 => n. *)

(* preclevel[m,n] chooses precedence levels based on the version of
   Mathematica in use: version < 3.0 => m; version >= 3.0 => n. *)

preclevel[m_,n_] := 
  If[$VersionNumber < 3.0, m, n];

(* TimesFormat implements the formatting of products *)

HoldPattern[TimesFormat[{x___},form___]] := 
  Module[{factors,constants,tensors,others},
  factors    = holdlist[x];
  (* constants = factors with no tensors in them *)
  constants  = Select[factors, ConstantQ];
  (* tensors = factors with rank > 0 *)
  tensors    = Select[factors, rank[#] > 0 & ];
  (* others = everything else *)
  others     = Select[factors, !ConstantQ[#] && !TrueQ[rank[#]>0] & ];
  (* unwrap the HoldForm from numbers and integer powers, so Mathematica
     can do its usual formatting of them. *)
  constants  = constants //. {HoldForm[y_^p_Integer] :> HoldForm[y]^p,
                              HoldForm[i_Integer] :> i,
                              HoldForm[i_Rational] :> i,
                              HoldForm[i_Real] :> i,
                              HoldForm[i_Complex] :> i};
  (* Now print the result: first constants, then others, then tensors. *)
  Which[
    constants==={},     TensorTimesFormat[others,tensors,form],
    constants==={-1},
                        Infix[{"-",
                            TensorTimesFormat[others,tensors,form]},"",
                            preclevel[141,311],None ],
    {others,tensors}==={{},{}},
                        Times @@ constants,
    True,               TensorTimesFormat[Join[{Times@@constants},others],tensors,form]
    ]
  ];

(************ Subroutines for TimesFormat *******************************)

(* holdlist [x,y,...] creates a list of the form {HoldForm[x],HoldForm[y],...}
   without evaluating any of its arguments. *)

holdlist[x_,y___] := Join[{HoldForm[x]},holdlist[y]];
holdlist[] := {}; 

(* TensorTimesFormat prints the "constant" and "other" factors and then
   calls SymmetricProductFormat to print the "tensor" factors. *)

HoldPattern[TensorTimesFormat[{},{s___},form___]] := 
  SymmetricProductFormat[{s},form];
HoldPattern[TensorTimesFormat[{c_},{},form___]] := Format[HoldForm[c],form];
HoldPattern[TensorTimesFormat[{c__},{s__},form___]] := 
  Format[ 
    Infix[ {c,SymmetricProductFormat[{s},form]}," ",preclevel[150,400],None],
    form];
HoldPattern[TensorTimesFormat[{c__},{},form___]] := 
  Format[ Infix[ {c}," ", preclevel[150,400],None], form ];

(* SymmetricProductFormat prints tensor factors with rank > 0 with 
   "*"'s in between *)

HoldPattern[SymmetricProductFormat[{x_},form___]] := 
  Format[HoldForm[x],form];
HoldPattern[SymmetricProductFormat[{x_,y__},form___]] := 
  Format[
    Infix[{x,y}, " * ", preclevel[152,402],None],
    form];
HoldPattern[SymmetricProductFormat[{},___]] := "";

(******************* TensorPowerFormat *****************************)

(* This prints tensor^p in exponential notation, instead of converting
   to a fraction *)

HoldPattern[TensorPowerFormat[Tensor[t_,{},{}],p_,form___]] := 
  Infix[{Format[HoldForm[Tensor[t,{},{}]],form],
         Superscript[Format[HoldForm[p],form]]},"",preclevel[190,590],NonAssociative];
HoldPattern[TensorPowerFormat[Tensor[t_,{i___},{j___}],p_,form___]] := 
  Infix[{SequenceForm["(",
                      Format[HoldForm[Tensor[t,{i},{j}]],form],
                      ")"],
         Superscript[Format[HoldForm[p],form]]},"",
         preclevel[190,590],NonAssociative];
HoldPattern[TensorPowerFormat[x_,p_,form___]] :=
  Infix[{Format[HoldForm[x],form],
         Superscript[Format[HoldForm[p],form]]},"",
         preclevel[190,590],NonAssociative];

(************************** TensorProduct **************************)

TensorProduct/: HoldPattern[Format[TensorProduct[x_,y_,z__]]] := 
                Infix[ {HoldForm[x],HoldForm[TensorProduct[y,z]]}, " (X) ", 
                preclevel[151,401]] /;
                $TensorFormatting;

TensorProduct/: HoldPattern[Format[TensorProduct[x_,y_]]] := 
                Infix[ {HoldForm[x],HoldForm[y]}, " (X) ", 
                preclevel[151,401]] /;
                $TensorFormatting;

TensorProduct/: HoldPattern[Format[TensorProduct[x_,y_,z__],TeXForm]] := 
                Infix[ {Format[HoldForm[x],TeXForm],
                        Format[HoldForm[TensorProduct[y,z]],TeXForm]}, 
                        "\\otimes ", preclevel[151,401]] /;
                $TensorTeXFormatting;

TensorProduct/: HoldPattern[Format[TensorProduct[x_,y_],TeXForm]] := 
                Infix[ {Format[HoldForm[x],TeXForm],
                        Format[HoldForm[y],TeXForm]}, 
                        "\\otimes ", preclevel[151,401]] /;
                $TensorTeXFormatting;

(************************** Wedge **************************)

Wedge/: HoldPattern[Format[Wedge[x_,y_,z__]]] := 
  Infix[{HoldForm[x],HoldForm[Wedge[y,z]]}, " ^ ", 
    preclevel[152,402]] /;
  $TensorFormatting;

Wedge/: HoldPattern[Format[Wedge[x_,y_]]] := 
  Infix[{HoldForm[x],HoldForm[y]}, " ^ ", 
    preclevel[152,402]] /;
  $TensorFormatting;

Wedge/: HoldPattern[Format[Wedge[x_,y_,z__],TeXForm]] := 
  Infix[{Format[HoldForm[x],TeXForm],
         Format[HoldForm[Wedge[y,z]],TeXForm]}, 
         "\\wedge ", preclevel[152,402]] /;
  $TensorTeXFormatting;

Wedge/: HoldPattern[Format[Wedge[x_,y_],TeXForm]] := 
  Infix[{Format[HoldForm[x],TeXForm],
         Format[HoldForm[y],TeXForm]}, 
         "\\wedge ", preclevel[152,402]] /;
  $TensorTeXFormatting;

(********************* Del  *******************************)

Del/: Format[Del,TeXForm] = "\\Del ";

Del/: HoldPattern[Format[Del[v_, x_, opt_ -> con_]]] := 
      SequenceForm[Format[Del],
                   ColumnForm[{
                     ColumnForm[{
                       SequenceForm["(",HoldForm[con],")"],
                       " "}, Left, Above],
                     HoldForm[v]
                     }, Left, Below],
                   " [",HoldForm[x],"]"] /;
      $TensorFormatting;

Del/: HoldPattern[Format[Del[v_, x_, opt_ -> con_],TeXForm]] := 
      SequenceForm[Format[Del,TeXForm],
                   Subscript[Format[HoldForm[v],TeXForm]],
                   Superscript[
                     SequenceForm["(",Format[HoldForm[con],TeXForm],")"]
                     ],
                   "\\left(",Format[HoldForm[x],TeXForm],"\\right)"] /;
      $TensorTeXFormatting;

Del/: HoldPattern[Format[Del[x_, opt_ -> con_]]] := 
      SequenceForm[Format[Del],
                   Superscript[SequenceForm["(",HoldForm[con],")"]],
                   " [",HoldForm[x],"]"] /;
      $TensorFormatting;

Del/: HoldPattern[Format[Del[x_, opt_ -> con_],TeXForm]] := 
      SequenceForm[Format[Del,TeXForm],
                   Superscript[
                     SequenceForm["(",Format[HoldForm[con],TeXForm],")"]
                     ],
                   "\\left(",Format[HoldForm[x],TeXForm],"\\right)"] /;
      $TensorTeXFormatting;

Del/: HoldPattern[Format[ Del[ Tensor[Basis,{L[i_]},{}], x_ ] ]] :=
      SequenceForm[ Format[Del],
                    IndexFormat[{L[i]},{}],
                    "[", HoldForm[x], "]" ] /;
      $TensorFormatting;

Del/: HoldPattern[Format[ Del[ Tensor[Basis,{L[i_]},{}], x_ ], TeXForm ]] :=
      SequenceForm[ TensorTeXFormat[ Del,
                                     {L[i]},{} ],
                    "\\left(",Format[HoldForm[x],TeXForm],"\\right)"] /;
      $TensorTeXFormatting;

Del/: HoldPattern[Format[Del[v_,x_]]] := 
      SequenceForm[Format[Del],
                   Subscript[HoldForm[v]],
                   " [",HoldForm[x],"]"] /;
      $TensorFormatting;

Del/: HoldPattern[Format[Del[v_,x_],TeXForm]] := 
      SequenceForm[Format[Del,TeXForm],
                   Subscript[Format[HoldForm[v],TeXForm]],
                   "\\left(",Format[HoldForm[x],TeXForm],"\\right)"] /;
      $TensorTeXFormatting;

(********************* Grad  *******************************)

Grad/: Format[Grad,TeXForm] = "\\grad ";

Grad/: HoldPattern[Format[Grad[x_, opt_ -> con_]]] := 
      SequenceForm[Format[Grad],
                   Superscript[SequenceForm["(",HoldForm[con],")"]],
                   " [",HoldForm[x],"]"] /;
      $TensorFormatting;

Grad/: HoldPattern[Format[Grad[x_, opt_ -> con_],TeXForm]] := 
      SequenceForm[Format[Grad,TeXForm],
                   Superscript[
                     SequenceForm["(",Format[HoldForm[con],TeXForm],")"]
                     ],
                   "\\left(",Format[HoldForm[x],TeXForm],"\\right)"] /;
      $TensorTeXFormatting;

(********************* Lie  *******************************)

Lie/: HoldPattern[Format[Lie[v_,x_]]] := 
      SequenceForm[Format[Lie],Subscript[HoldForm[v]]," [",HoldForm[x],"]"] /;
      $TensorFormatting;

Lie/: HoldPattern[Format[Lie[v_,x_],TeXForm]] := 
      SequenceForm[Format[Lie,TeXForm],
                   Subscript[Format[HoldForm[v],TeXForm]],
                   "\\left(",Format[HoldForm[x],TeXForm],"\\right)"] /;
      $TensorTeXFormatting;

Lie/: Format[Lie,TeXForm] = "\\Lie ";

(********************* Summation  *******************************)

Summation/: HoldPattern[Format[Summation[x___]]] :=
                    SequenceForm["(",Format[HoldForm[x]],")"] /;
                    $TensorFormatting;

Summation/: HoldPattern[Format[Summation[x___],TeXForm]] :=
                    SequenceForm["\\left(",Format[HoldForm[x],TeXForm],
                                 "\\right)"] /;
                    $TensorTeXFormatting;

(********************* Div  *******************************)

Div/: Format[Div,TeXForm] = "\\div ";

Div/: HoldPattern[Format[Div[t_,opt_->met_]]] :=
  SequenceForm[Format[Div],
               Superscript["("],
               Superscript[HoldForm[met]],
               Superscript[")"],
               "[",HoldForm[t],"]"] /;
  $TensorFormatting;

Div/: HoldPattern[Format[Div[t_,opt_->met_],TeXForm]] :=
  SequenceForm[Format[Div,TeXForm],
               Superscript["("],
               Superscript[Format[HoldForm[met],TeXForm]],
               Superscript[")"],
               "\\left(",Format[HoldForm[t],TeXForm],"\\right)"] /;
  $TensorTeXFormatting;

(********************* Alt, Sym  *******************************)

Alt/: Format[Alt,TeXForm] = "\\Alt ";

Sym/: Format[Sym,TeXForm] = "\\Sym ";

(********************* Int  *******************************)

Int/: HoldPattern[Format[Int[x_,y_]]] :=
      SequenceForm[Format[Int],
                   Subscript[HoldForm[x]],
                   " [",HoldForm[y],"]"] /;
      $TensorFormatting;

Int/: HoldPattern[Format[Int[x_,y_],TeXForm]] :=
      SequenceForm[Format[Int,TeXForm],
                   Subscript[Format[HoldForm[x],TeXForm]],
                   "\\left(",Format[HoldForm[y],TeXForm],"\\right)"] /;
      $TensorTeXFormatting;

Int/: Format[Int,TeXForm] = "i";

(********************** Inner ***************************************)

Unprotect[Inner];

Inner/: HoldPattern[Format[Inner[x_,y_]]] := 
  SequenceForm["<",HoldForm[x],", ",HoldForm[y],">"];

Inner/: HoldPattern[Format[Inner[x_,y_],TeXForm]] := 
  SequenceForm["\\langle ",HoldForm[x],",",HoldForm[y],"\\rangle "];

Protect[Inner];

(********************** Inverse **********************************)

Unprotect[Inverse];

Inverse/: HoldPattern[Format[Inverse[Tensor[x_,{},{}]]]] := 
          Infix[ {HoldForm[x], Superscript[-1]}, "", preclevel[180,580] ];

Inverse/: HoldPattern[Format[Inverse[x_]]] :=
          Infix[ {SequenceForm[ "(",HoldForm[x],")" ],
                  Superscript[-1]}, "", preclevel[180,580] ];

Inverse/: HoldPattern[Format[Inverse[Tensor[x_,{},{}]],TeXForm]] := 
          Infix[ {Format[HoldForm[x],TeXForm], Superscript[-1]}, "", 
                 preclevel[180,580] ];

Inverse/: HoldPattern[Format[Inverse[x_],TeXForm]] :=
          Infix[ {SequenceForm[ "(",Format[HoldForm[x],TeXForm],")"],
                  Superscript[-1]}, "", preclevel[180,580] ];

Protect[Inverse];

(********************** HodgeInner **********************************)

HodgeInner/: HoldPattern[Format[HodgeInner[x_,y_]]] := 
  SequenceForm["<<",HoldForm[x],", ",HoldForm[y],">>"];

HodgeInner/: HoldPattern[Format[HodgeInner[x_,y_],TeXForm]] := 
  SequenceForm["\\langle\\langle ",HoldForm[x],",",
                                   HoldForm[y],"\\rangle\\rangle "];

(********************* Power *******************************)

Unprotect[Power];

Power/: HoldPattern[Format[t_^p_Integer]] := 
  TensorPowerFormat[t,p] /; 
  !FreeQ[Hold[t],Tensor[_,_,_]] && $TensorFormatting;
Power/: HoldPattern[Format[t_^p_Integer,TeXForm]] := 
  TensorPowerFormat[t,p,TeXForm] /; 
  !FreeQ[Hold[t],Tensor[_,_,_]] && $TensorTeXFormatting;

Protect[Power];

(********************* Extd, ExtdStar  *******************************)

Extd/: Format[Extd] := "d" /; $TensorFormatting;

Extd/: Format[Extd,TeXForm] := "d" /; $TensorTeXFormatting;

ExtdStar/: Format[HoldPattern[ExtdStar[t_]]] := 
  SequenceForm["d",Superscript["*"],"[",HoldForm[t],"]"] /;
  $TensorFormatting;

ExtdStar/: Format[HoldPattern[ExtdStar[t_]],TeXForm] := 
  SequenceForm["d^*\\left(",Format[t,TeXForm],"\\right)"] /;
  $TensorTeXFormatting;

(********************* Conjugate *******************************)

Unprotect[Conjugate];

Conjugate/: HoldPattern[Format[Conjugate[x_Symbol]]] := 
  BarFormat[HoldForm[x]];

Conjugate/: HoldPattern[Format[Conjugate[x_Symbol], TeXForm]] := 
  BarFormat[HoldForm[x],TeXForm];

Protect[Conjugate];

(********************* ERROR  *******************************)

ERROR/: Format[x_ERROR] := InputForm[x];

(*************************** Index *****************************)

(* Functions for determining index properties *)

IndexName[L[i_]] := i;
IndexName[U[i_]] := i;

IndexAltitude[L[_]] := Lower;
IndexAltitude[U[_]] := Upper;

IndexQ[] := True;
IndexQ[L[_]] := True;
IndexQ[U[_]] := True;
IndexQ[_] := False;
IndexQ[i_,j__] := And @@ (IndexQ/@{i,j});

LowerIndexQ[L[_]]    := True;
UpperIndexQ[U[_]]    := True;
LowerIndexQ[_] := False;
UpperIndexQ[_] := False;

Pair[L[i_]] := U[i];
Pair[U[i_]] := L[i];

PairQ[i_,j_] := (Pair[i]===j);

Lower[L[i_]] := L[i];
Lower[U[i_]] := Conjugate[L[i]];
Lower[x_]    := x;
Raise[L[i_]] := Conjugate[U[i]];
Raise[U[i_]] := U[i];
Raise[x_]    := x;

(* Conjugation rules for indices *)

L/: Conjugate[L[x_[Bar]]] := L[x];
L/: Conjugate[L[x_]]      := L[x] /; Type[Bundle[x]] === Real;
L/: Conjugate[L[x_]]      := L[x[Bar]] ;
U/: Conjugate[U[x_[Bar]]] := U[x];
U/: Conjugate[U[x_]]      := U[x] /; Type[Bundle[x]] === Real;
U/: Conjugate[U[x_]]      := U[x[Bar]] ;

(* Bundle[indexname] returns the name of the bundle associated with
   indexname.  For index names defined by the user, it is set by
   DefineIndex.  For dummy index names such as i1, i2 etc. generated
   by Ricci, it is set when the index name is created. *)

Bundle[x_[Bar]] := Conjugate[Bundle[x]];
Bundle[L[i_]] := Bundle[i];
Bundle[U[i_]] := Bundle[i];

(* If the user uses an index name like i37 which has never been used
   before, associate it with the same bundle as i. *)

Bundle[i_Symbol] := Module[{prefix},
  prefix = ToExpression[StringJoin @@
           (Select[Characters[ToString[i]],LetterQ]) ];
  If[prefix===i || Bundle[prefix]===Null,
     (*then*) Return[Null],
     (*else*) i/: Bundle[i] = Bundle[prefix];
              Protect[i];
              Return[Bundle[prefix]]
    ]
  ];

Bundle[_] := Null;

(* LB & UB are abbreviations for conjugate indices. *)

LB[i_] := Conjugate[L[i]];
UB[i_] := Conjugate[U[i]];

(************************ Index ordering functions ***********************)

(* We order indices in our own way, instead of using Mathematica's ordering.
   First by name, then by altitude. *)

IndexOrderedQ[ {a_[b_],c_[d_]} ] := OrderedQ[ {{b,a},{d,c}} ];
IndexOrderedQ[{i:(L[_]|U[_])..}] := OrderedQ[ {#[[1]],#}& /@ {i} ];
IndexOrderedQ[{x___}] := OrderedQ[{x}];

IndexOrderedQ[i___] := IndexOrderedQ[{i}];

(* The following form is used by DefineTensor (actually by the subroutine 
   defTenSym) to define arbitrary tensor symmetries.  It returns True if the
   list i is better-ordered than the list j, according to our rule for 
   ordering indices. *)

IndexOrderedQ[{i__},{j__}] :=
  OrderedQ[ { {#[[1]],#}& /@ {i}, {#[[1]],#}& /@ {j} } ];

IndexSort[ {i:(L[_]|U[_])..} ] := Last /@ Sort[ {#[[1]],#}& /@ {i} ];
IndexSort[ {x___} ] := Sort[{x}];

(* This function is analogous to Mathematica's Signature, except instead of 
   Mathematica's canonical order it uses our rule for ordering indices. *)

IndexSignature[list_] := 
  Signature[ Position[ IndexSort[list], #]& /@ list];

(******************** GetIndices, etc. *************************************)

GetIndices[exp_] := Occurrences[ exp, L[_Symbol]|U[_Symbol]|
                                      L[_Symbol[Bar]]|U[_Symbol[Bar]] ];

PairedQ[x_,y_] := !Intersection[GetIndices[x], Pair /@ GetIndices[y]] === {};

GetFreeIndices[a_Plus] := 
  GetFreeIndices[ First[a] ];

GetFreeIndices[a_] := 
  GetIndices[a] //. 
    { {x___,L[i_],y___,U[i_],z___} :> {x,y,z} /; Dimension[Bundle[i]]=!=1,
      {x___,U[i_],y___,L[i_],z___} :> {x,y,z} /; Dimension[Bundle[i]]=!=1 };

(* GetDummyIndices returns the lower indices of dummy pairs.  Indices
   associated with one-dimensional bundles are all considered free. *)

GetDummyIndices[a_] := 
  Select[ 
    Intersection[
           Occurrences[ a, L[_Symbol]|L[_Symbol[Bar]] ],
           Pair /@ Occurrences[ a, U[_Symbol]|U[_Symbol[Bar]] ] ],
    Dimension[Bundle[#]] =!= 1 & ];

(************************ DefineIndex ****************************)

Attributes[DefineIndex] = {Listable};

Options[DefineIndex] := {Quiet -> $Quiet, TeXFormat -> Null};

DefineIndex[]  := Message[DefineIndex::argm,"DefineIndex",0,2];
DefineIndex[_] := Message[DefineIndex::argmu,"DefineIndex",2];

DefineIndex[i_, bundle_, opts___] := Module[{quiet,tex,protectedbun},
  If[ !optionsOK[DefineIndex,opts], Return[] ];
  {quiet,tex} = {Quiet,TeXFormat} /. {opts} /. Options[DefineIndex];
  Which[
    !newsymQ[i] && !(IndexQ[i] && Bundle[i]===bundle),
      Message[DefineIndex::invalid, i, "index name"],
    !BundleQ[bundle] && !newsymQ[bundle],
      Message[DefineIndex::invalid, bundle, "bundle name"],
    Head[tex] =!= String && tex =!= Null,
      Message[DefineIndex::invalid, tex, "TeX format"],
    True,
      Unprotect[i];
      i/: IndexQ[i] = True;
      i/: Bundle[i] = bundle;
      protectedbun = Unprotect[bundle];
      bundle/: BundleIndices[bundle] = 
                          Union[ BundleIndices[bundle], {i} ];
      Protect[Evaluate[protectedbun]];
      Protect[i];
      DeclareIndex[i, TeXFormat -> tex];
      If[ !quiet, (* then print message *)
        Print["Index ",i," associated with ",bundle];
        ];
      Return[i];
    ];
  ];

(************************* DeclareIndex *****************************)

(* This function is called by Declare when the first argument is
   an index. See Constant.m *)

Options[DeclareIndex]  := {TeXFormat -> Null};

DeclareIndex[ i_, opts___ ] := Module[{bun,tex},
  {tex} = ({TeXFormat} /. {opts} /. Options[DeclareIndex]);
  Which[
    Head[tex]=!=String && tex =!= Null,
      Message[DefineIndex::invalid, tex, "TeXFormat"],
    True,
      Unprotect[i];
      If[ tex =!= Null,
          i/: Format[i,TeXForm] = tex ];
      Protect[i];
      Return[i];
    ]
  ];
(************************* UndefineIndex *****************************)

Attributes[UndefineIndex] = {Listable};

Options[UndefineIndex] := {Quiet -> $Quiet};

UndefineIndex[i_,opts___] := Module[{quiet,bun},
  If[!optionsOK[UndefineIndex,opts], Return[]];
  quiet = Quiet /. {opts} /. Options[UndefineIndex];
  bun = Bundle[i];
  If[bun===Null, 
     (*then error *) 
        Print["ERROR: ",i," is not an index"];
        Return[],
     (*else OK*)     
        Unprotect[i];
        ClearAll[i];
        Unprotect[Release[bun]];
        Release[bun]/: BundleIndices[bun] = Complement[BundleIndices[bun],{i}];
        Protect[Release[bun]];
        If[!quiet, 
          Print["Index ",i," is no longer associated with ",bun,"."]]
    ];
  ];

(***************** Functions for generating new dummy indices **********)

(* NewBundleSymbol[bundle] just generates a new index name, without
   L or U wrapped around it. *)

NewBundleSymbol[Null] := Null;

NewBundleSymbol[Conjugate[bundle_]] := NewBundleSymbol[bundle][Bar];

NewBundleSymbol[bundle_] := BundleDummy[bundle] /; Dimension[bundle]===1;

NewBundleSymbol[bundle_] := Module[{i},
  i = Unique[ToString[BundleDummy[bundle]]];
  Bundle[i] ^= bundle;
  Protect[i];
  Return[i];
  ];

(* NewLowerIndex & NewUpperIndex generate new indices of the form
   L[i1] or U[i1]. *)

NewLowerIndex[Conjugate[bundle_]] := L[NewBundleSymbol[bundle][Bar]];
NewUpperIndex[Conjugate[bundle_]] := U[NewBundleSymbol[bundle][Bar]];
NewLowerIndex[bundle_] := L[NewBundleSymbol[bundle]];
NewUpperIndex[bundle_] := U[NewBundleSymbol[bundle]];

(* NewIndex[ {listofbundles, variance} ] generates a list of upper 
   (Contravariant) or lower (Covariant) index names corresponding
   to the bundles in "listofbundles". *)

NewIndex[{{Any},variance_}] := 
  NewIndex[{Flatten[{$DefaultTangentBundle}],variance}];

NewIndex[{{bundles__},variance_}] :=
  NewIndex[#,variance]& /@ {bundles};

NewIndex[bundle_,variance_] := 
  If[variance===Covariant,
     (*then*) NewLowerIndex[bundle],
     (*else*) NewUpperIndex[bundle] ];

(*************************** Products *****************************)

(******************* TensorProduct ***************************************)

TensorProduct[x_] := x;

(* Associativity is implemented by the following rule.  We do not use
   the Flat attribute, because there seems to be a bug in Mathematica (1.2)
   applying certain rules to Flat functions. *)

HoldPattern[TensorProduct[w___, TensorProduct[x__], z___]] := 
  TensorProduct[w,x,z];

(* Multilinearity. *)

HoldPattern[TensorProduct[ w___, (x_ /; rank[x]===0) * (y_), z___ ]] := 
  x * TensorProduct[w,y,z];

HoldPattern[TensorProduct[ w___, (x_) + (y_), z___ ]] := 
  TensorProduct[w,x,z] + TensorProduct[w,y,z];

(* Factor out scalars *)

HoldPattern[TensorProduct[w___, (x_ /; rank[x]===0) , z___]] := 
  x * TensorProduct[w,z];

TProd = TensorProduct;

TensorProduct/: Conjugate[HoldPattern[TensorProduct[x__]]] := 
                  Conjugate /@ TensorProduct[x];

(* Here's where we insert indices *)

TensorProduct/: TCompute[HoldPattern[TensorProduct[ x_, y__ ]], {i___}] := 
  Return[ InsertIndices[ x, Take[{i},rank[x]] ] * 
          InsertIndices[ TensorProduct[y], Drop[{i},rank[x]] ] 
        ];

(*********************** SymmetricProduct **********************)

(* SymmetricProduct is represented internally by ordinary 
   multiplication. *)

HoldPattern[SymmetricProduct[x___]] := Times[x];

(*********************** Times *********************************)

(* Insertion of indices into products is done here. *)

Unprotect[Times];

Times/: TCompute[ x_Times, {} ] := InsertIndices[#,{}]& /@ x;

Times/: TCompute[ HoldPattern[Times[x__]], {i__} ] := 
  If[ !(And @@ SymmetricQ /@ {x}),
      (*then*) Print["ERROR: non-symmetric tensor used in symmetric product"];
               ERROR[SymmetricProduct[x][i]],
      (*else*) TCompute[ Unevaluated[ Sym[ TensorProduct[x] ] ], {i} ]
               (* Unevaluated is here so that Sym[ a(X)b ] won't get
                  transformed back to a*b, causing an infinite loop *)
    ];

Protect[Times];

(****************** Inner *********************************)

Unprotect[Inner];
Off[Inner::argb];

(* Bilinearity *)
Inner[x_ + y_, z_] := Inner[x,z] + Inner[y,z];
Inner[(c_ /; rank[c]===0) * x_, z_] := c Inner[x,z];
Inner[x_, y_ + z_] := Inner[x,y] + Inner[x,z];
Inner[x_, (c_ /; rank[c]===0) * y_] := c Inner[x,y];

Inner[_, 0] := 0;
Inner[0, _] := 0;

Inner[x_, y_] := (
  Print["ERROR: ",x," and ",y," don't have the same rank"];
  ERROR[Inner[x,y]]
  ) /; !Rank[x]===Rank[y];

Inner[Tensor[Basis,{i_},{}], x_] := InsertIndices[x, {i}];
Inner[x_, Tensor[Basis,{i_},{}]] := InsertIndices[x, {i}];

Inner[x_, x_] := Plus @@ Dimension /@ Bundles[x][[1]] /; MetricQ[x];

Inner[x_, y_] := x y /; rank[x]===0;

Inner[x_, y_] := 0 /; rank[x] > 1 && 
                         ((SymmetricQ[x] && SkewQ[y]) ||
                             (SkewQ[x] && SymmetricQ[y]));

(* If Inner represents a pure pairing between covariant and contravariant
   tensors, make sure the covariant one comes first. *)

Inner[x_, y_] := Inner[y, x] /; Union[Variance[y]]==={Covariant} &&
                                Union[Variance[x]]==={Contravariant};

(* Otherwise, put them in canonical order. *)

Inner[x_, y_] := Inner[y, x] /; !OrderedQ[{x,y}] && 
                                (!Union[Variance[x]]==={Covariant} ||
                                 !Union[Variance[y]]==={Contravariant});

Inner[ x___, HoldPattern[Grad[y_]], z___ ] := Inner[x,Del[y],z];

Inner[ HoldPattern[Extd[u_]], v_ ] := Del[v,u] /; VectorFieldQ[v];

Inner/: Conjugate[HoldPattern[Inner[x_,y_]]] :=
                  Conjugate /@ Inner[x,y];

Inner/: TCompute[ HoldPattern[Inner[x_, y_]], {} ] := (
    Plus @@ Flatten [ Outer @@ Join[
      {InsertIndices[x, Pair /@ {##}] * InsertIndices[y, {##}] &},
      NewIndex /@ Transpose[{Bundles[y],Variance[y]}] ]]
    );

Protect[Inner];

(***************** Norm *********************************)

Unprotect[Norm];
Norm[x_] := Inner[x,Conjugate[NewDummy[x]]]^(1/2);
Protect[Norm];

(****************** HodgeInner *********************************)

(* bilinearity *)
HodgeInner[x_ + y_, z_] := HodgeInner[x,z] + HodgeInner[y,z];
HodgeInner[(c_ /; rank[c]===0) * x_, z_] := c HodgeInner[x,z];
HodgeInner[x_, y_ + z_] := HodgeInner[x,y] + HodgeInner[x,z];
HodgeInner[x_, (c_ /; rank[c]===0) * y_] := c HodgeInner[x,y];

HodgeInner[_, 0] := 0;
HodgeInner[0, _] := 0;

HodgeInner[x_, y_] := (
  Print["ERROR: ",x," and ",y," don't have the same rank."];
  ERROR[HodgeInner[x,y]]
  ) /; !rank[x]===rank[y];

HodgeInner[x_, y_] := (
  Print["ERROR: HodgeInner applied to a non-skew tensor."];
  ERROR[HodgeInner[x,y]]
  ) /; !SkewQ[x] || !SkewQ[y];

HodgeInner[x_, y_] := x y /; rank[x]===0;

HodgeInner[x_, y_] := Inner[x,y] /; rank[x]===1;

(* If HodgeInner represents a pure pairing between covariant and contravariant
   tensors, make sure the covariant one comes first. *)

HodgeInner[x_, y_] := HodgeInner[y, x] /; Union[Variance[y]]==={Covariant} &&
                                Union[Variance[x]]==={Contravariant};

(* Otherwise, put them in canonical order. *)

HodgeInner[x_, y_] := HodgeInner[y, x] /; !OrderedQ[{x,y}] && 
                                (!Union[Variance[x]]==={Covariant} ||
                                 !Union[Variance[y]]==={Contravariant});

HodgeInner/: Conjugate[HoldPattern[HodgeInner[x_,y_]]] :=
                  Conjugate /@ HodgeInner[x,y];

(* Here's where we insert indices. *)

HodgeInner/: TCompute[ HoldPattern[HodgeInner[x_, y_]], {} ] := 
  IntFactor[rank[x],rank[y]] * InsertIndices[ Inner[x,y], {} ];

(***************** HodgeNorm *********************************)

HodgeNorm[x_] := HodgeInner[x,Conjugate[NewDummy[x]]]^(1/2);

(****************** Int *********************************)

(* IntFactor is a numerical factor, depending on WedgeConvention,
   which is applied when Int is computed. *)

IntFactor[m_,n_] := 
   n!/((n-m)! * WedgeFactor[m]^2 * WedgeFactor[{m,n-m}]);

(* Bilinearity. *)

Int[x_ + y_, z_] := Int[x,z] + Int[y,z];
Int[(c_ /; rank[c]===0) * x_, z_] := c Int[x,z];
Int[x_, y_ + z_] := Int[x,y] + Int[x,z];
Int[x_, (c_ /; rank[c]===0) * y_] := c Int[x,y];

Int[x_, y_] := (
  Print["ERROR: Int applied to a non-skew tensor."];
  ERROR[Int[x,y]]
  ) /; !SkewQ[x] || !SkewQ[y];

(* Antiderivation rule. *)

HoldPattern[Int[(v_ /; rank[v]===1), Wedge[y_,z__]]] :=
  Wedge[ Int[v,y], z ] + (-1)^rank[y] Wedge[ y, Int[v,Wedge[z]] ];

(* v _| v _| x = 0  if v is a vector field and x is a form. *)

HoldPattern[Int[(v_ /; rank[v]===1), Int[v_, x_] ] ] := 0;

Int[x_, y_] := 0   /; rank[x] > rank[y];
Int[x_, y_] := x y /; rank[x]===0;

Int[x_, y_] := HodgeInner[x,y] /; rank[x]===rank[y];

Int/: Conjugate[HoldPattern[Int[x_,y_]]] :=
                  Conjugate /@ Int[x,y];

(* Here's where we insert indices. *)

Int/: TCompute[ HoldPattern[Int[x_, y_]], {i___} ] := (
  Print["ERROR: non-skew tensor used in Int"];
  ERROR[Int[x,y][i]]
  ) /; (!SkewQ[x] || !SkewQ[y]);

Int/: TCompute[ HoldPattern[Int[x_, y_]], {i___} ] := (
  IntFactor[rank[x],rank[y]] *
    Plus @@ Flatten [ Outer @@ Join[
      {InsertIndices[x, {##}] * InsertIndices[y, Join[ Pair /@ {##}, {i} ] ] &},
      NewIndex /@ Transpose[{Bundles[x],Variance[x]}] ]]
  );

(******************** Dot *******************************************)

Unprotect[Dot];

Dot[w___, x_ + y_, z___] := Dot[w,x,z] + Dot[w,y,z];
Dot[w___, (f_ /; rank[f]===0) * x_, z___] := f Dot[w,x,z];
Dot[___, 0, ___] := 0;

Dot[ a___, HoldPattern[TensorProduct[x_,y__]], b___ ] := 
  TensorProduct[ Dot[a,x], Dot[TensorProduct[y],b] ];

Dot[ a___, HoldPattern[Grad[b_]], c___ ] := Dot[a,Del[b],c] /;
  Length[{c}] =!= 0 || Rank[b] === 0;

Dot[ a___, HoldPattern[Del[f_,con___Rule]], (v_ /; Rank[v]===1)] := 
  Dot[ a, Del[v,f,con] ];

Dot[x_, y_] := Inner[x,y] /; rank[x]===rank[y]===1 &&
                             Head[x] =!= Dot &&
                             Head[y] =!= Dot;

(* Kronecker acts as identity in matrix multiplication. *)

Dot[ a___, HoldPattern[Tensor[Kronecker,{},{}]], b___ ] := Dot[a,b] /;
  Length[{a,b}] >= 1;

Dot[ x_, y__, z_ ] := (
  Print["ERROR: Rank 1 tensor used incorrectly in Dot product."];
  ERROR[Dot[x,y,z]]
  ) /; MemberQ[ rank /@ {y}, 1 ] && !MemberQ[ Head /@ {x}, List ];

Dot[ x___ ] := (
  Print["ERROR: Rank 0 tensor used in Dot product."];
  ERROR[Dot[x]]
  ) /; MemberQ[ rank /@ {x}, 0 ] && !MemberQ[ Head /@ {x}, List ];

Dot/: Conjugate[HoldPattern[x_Dot]] := Conjugate /@ x;

Dot/: TCompute[ HoldPattern[Dot[x__]], {i___} ] := 
  Module[{bundles,vars,n,startpos},
  bundles = (Last /@ Bundles /@ Drop[{x},-1]) /. 
    {{Any} -> Flatten[{$DefaultTangentBundle}]};
  vars    = Last /@ Variance /@ Drop[{x},-1];
  startpos[1] = 1;
  startpos[2] = Rank[ First[{x}] ];
  startpos[n_] := startpos[n-1] + Rank[ {x}[[n-1]] ] - 2;
  Plus @@ Flatten [ Outer @@ Join[
    {InsertIndices[ First[{x}], 
                    Append[ Take[{i},Rank[First[{x}]]-1], #1]
                  ] *
     Product[ 
       InsertIndices[ {x}[[n]], 
                      Join[ {Pair[ {##}[[n-1]] ]},
                            Take[ {i}, 
                                  {startpos[n], startpos[n] + Rank[ {x}[[n]] ] - 3}
                                ], 
                            { {##}[[n]] }
                          ]
                    ],
       {n, 2, Length[{x}]-1}
       ] *
     InsertIndices[ Last[{x}],
                    Prepend[ Take[{i},-(Rank[Last[{x}]]-1)], Pair @ Last[{##}] ]
                  ] &
    },
    NewIndex /@ Transpose[{bundles,vars}] ]
    ]
  ];

Protect[Dot];

(***************************** Wedge *********************************)

WedgeFactor[{p___}] := 
  Switch[$WedgeConvention,
    Alt,   1,
    Det,   wfactor[p],
    _,     Message[WedgeConvention::invalid, $WedgeConvention, 
                   "$WedgeConvention"];
           1
  ];

(** WedgeFactor[p] is an abbreviation for WedgeFactor[{1,...,1}] **)

WedgeFactor[p_] := WedgeFactor[Table[1,{p}]];

wfactor[p___] := (Plus[p])!/(Times @@ Factorial /@ {p});

HoldPattern[Wedge[x_]] := x;

(* Associativity. *)

HoldPattern[Wedge[w___, Wedge[x__], z___]] := Wedge[w,x,z];

(* Multilinearity. *)

HoldPattern[Wedge[ w___, x_ + y_, z___ ]] := 
  Wedge[w,x,z] + Wedge[w,y,z];

HoldPattern[Wedge[ w___, (x_ /; rank[x]===0) * y_, z___ ]] := 
  x*Wedge[w,y,z];

(* Factor out zero forms. *)

HoldPattern[Wedge[w___, (x_ /; rank[x]===0), z___]] := x Wedge[w,z];

(* Anticommutativity. *)

HoldPattern[Wedge[w___, x_, y_, z___]] := 
  (-1) ^ (rank[x] * rank[y]) Wedge[w,y,x,z] /; !OrderedQ[{x,y}];

HoldPattern[Wedge[w___, x_, y___, x_, z___]] := 0 /; OddQ[rank[x]];

Wedge/:     Conjugate[HoldPattern[Wedge[x__]]] := Conjugate /@ Wedge[x];

(** Here's what happens when you put indices in. **)

Wedge/: TCompute[ HoldPattern[Wedge[x_,y__]], {i___} ] := 
  If[ !(And @@ SkewQ /@ {x,y}),
      (*then*) Print["ERROR: Non-skew tensor used in wedge product"];
               ERROR[Wedge[x,y][i]],
      (*else*) WedgeFactor[Rank/@{x,y}] * 
                 TCompute[ Unevaluated[Alt[TensorProduct[x,y]]], {i} ]
    ];

(*************************** Riemann *****************************)

(******************* makeRiemannian **************************)

(* This is called by DefineBundle whenever the option 
   MetricType -> Riemannian is specified. *)

makeRiemannian[bun_,rm_,rc_,sc_,rmconv_,quiet_] := Module[{rmsign},
  Unprotect[bun];
  rmsign = If[rmconv === SecondUp,
              (*then*)   1,
              (*else*)   -1];
  bun/: RiemannConvention[bun] = rmconv;
  (* Define the Riemannian curvature tensors *)
  DefineTensor[ rm, 4, Symmetries -> RiemannSymmetries, Bundle->bun,
                Riemannian -> True ];
  DefineTensor[ rc, 2, Symmetries -> Symmetric, Bundle->bun,
                Riemannian -> True ];
  DefineTensor[ sc, 0, Bundle -> bun, Riemannian -> True ];
  bun/: RiemannTensor[bun] = rm;
  bun/: RicciTensor[bun] = rc;
  bun/: ScalarCurv[bun] = sc;
  defCurvatureRelations[rm,rc,sc,bun,rmsign];
  Protect[bun];
  ];

(* defCurvatureRelations defines all the standard transformations that
   will take place on this bundle: Curv -> Rm, Rm -> Rc, Rc -> Sc. *)

defCurvatureRelations[Tensor[rm_,_,_],Tensor[rc_,_,_],Tensor[sc_,_,_],
                      bun_,sgn_] := (
      Unprotect[rm,rc,Curv];
      rm/: HoldPattern[Tensor[rm,{a_[i_],j_,b_[i_],l_},{x___}]] :=
           -sgn * CovD[rc[j,l],{x}] /; !a===b;
      rm/: HoldPattern[Tensor[rm,{i_,a_[j_],b_[j_],l_},{x___}]] :=
            sgn * CovD[rc[i,l],{x}] /; !a===b;
      rm/: HoldPattern[Tensor[rm,{i_,a_[l_],k_,b_[l_]},{x___}]] :=
           -sgn * CovD[rc[i,k],{x}] /; !a===b;
      rc/: HoldPattern[Tensor[rc,{a_[i_],b_[i_]},{x___}]] :=
           CovD[sc,{x}] /; !a===b;
      Curv/: HoldPattern[Tensor[Curv, {i_,j_,k_,l_}, {x___}]] :=
           sgn * CovD[rm[i,j,k,l],{x}] /;
           Bundle[i]==Bundle[j]==bun;
      Protect[rm,rc,Curv]
      );

(******************* Rules for Bianchi identities ********************)

(* FirstBianchiRule implements the first Bianchi identity, turning sums of
   two Riemann tensors into single terms when possible. *)

FirstBianchiRule = {
  a_. * HoldPattern[Tensor[rm_,{i_,j_,k_,l_},{x___}]] + 
  a_. * HoldPattern[Tensor[rm_,{i_,l_,j_,k_},{x___}]]
    :> a * Tensor[rm,{i,k,j,l},{x}] /; 
    (Symmetries /. TensorData[rm]) === RiemannSymmetries,
  a_. * HoldPattern[Tensor[rm_,{i_,j_,k_,l_},{x___}]] +
  b_. * HoldPattern[Tensor[rm_,{i_,k_,j_,l_},{x___}]]
    :> a * Tensor[rm,{i,l,k,j},{x}] /; 
    (Symmetries /. TensorData[rm]) === RiemannSymmetries && a+b==0,
  a_. * HoldPattern[Tensor[rm_,{i_,j_,k_,l_},{x___}]] + 
  b_. * HoldPattern[Tensor[rm_,{i_,l_,k_,j_},{x___}]]
    :> a * Tensor[rm,{i,k,j,l},{x}] /;     
    (Symmetries /. TensorData[rm]) === RiemannSymmetries && a+b==0
  };

(* SecondBianchiRule implements the second Bianchi identity similarly. *)

SecondBianchiRule = {
  a_. * HoldPattern[Tensor[rm_,{p_,q_,i_,j_},{k_,x___}]] + 
  a_. * HoldPattern[Tensor[rm_,{p_,q_,j_,k_},{i_,x___}]]
    :> a * Tensor[rm,{p,q,i,k},{j,x}] /; RiemannQ[rm],
  a_. * HoldPattern[Tensor[rm_,{p_,q_,i_,j_},{k_,x___}]] +
  b_. * HoldPattern[Tensor[rm_,{p_,q_,i_,k_},{j_,x___}]]
    :> a * Tensor[rm,{p,q,k,j},{i,x}] /; RiemannQ[rm] && a+b==0,
  a_. * HoldPattern[Tensor[rm_,{p_,q_,j_,k_},{i_,x___}]] + 
  b_. * HoldPattern[Tensor[rm_,{p_,q_,i_,k_},{j_,x___}]]
    :> a * Tensor[rm,{p,q,j,i},{k,x}] /; RiemannQ[rm] && a+b==0,
  a_. * HoldPattern[Tensor[rm_,{p_,q_,j_,k_},{i_,x___}]] + 
  b_. * HoldPattern[Tensor[rm_,{i_,k_,p_,q_},{j_,x___}]]
    :> a * Tensor[rm,{p,q,j,i},{k,x}] /; RiemannQ[rm] && a+b==0,
  a_. * HoldPattern[Tensor[rm_,{p_,q_,i_,j_},{k_,x___}]] + 
  a_. * HoldPattern[Tensor[rm_,{j_,k_,p_,q_},{i_,x___}]]
    :> a * Tensor[rm,{p,q,i,k},{j,x}] /; RiemannQ[rm],
  a_. * HoldPattern[Tensor[rm_,{i_,j_,p_,q_},{k_,x___}]] + 
  a_. * HoldPattern[Tensor[rm_,{j_,k_,p_,q_},{i_,x___}]]
    :> a * Tensor[rm,{p,q,i,k},{j,x}] /; RiemannQ[rm],
  a_. * HoldPattern[Tensor[rm_,{i_,j_,p_,q_},{k_,x___}]] +
  b_. * HoldPattern[Tensor[rm_,{i_,k_,p_,q_},{j_,x___}]]
    :> a * Tensor[rm,{p,q,k,j},{i,x}] /; RiemannQ[rm] && a+b==0,
  a_. * HoldPattern[Tensor[rm_,{j_,k_,p_,q_},{i_,x___}]] + 
  b_. * HoldPattern[Tensor[rm_,{i_,k_,p_,q_},{j_,x___}]]
    :> a * Tensor[rm,{p,q,j,i},{k,x}] /; RiemannQ[rm] && a+b==0
  };

(* ContractedBianchiRules implements the once- and twice-contracted
   second Bianchi identity, always replacing Rm by Rc and Rc by their
   contractions (which will get transformed automatically to Rc and
   Sc) whenever possible. *)

ContractedBianchiRules = {
  HoldPattern[Tensor[rm_,{a_[i_],j_,k_,l_},{b_[i_],x___}]] :>
       Tensor[rm,{a[i],j,b[i],l},{k,x}] + 
       Tensor[rm,{a[i],j,k,b[i]},{l,x}] /; RiemannQ[rm],
  HoldPattern[Tensor[rm_,{i_,a_[j_],k_,l_},{b_[j_],x___}]] :>
       Tensor[rm,{i,a[j],b[j],l},{k,x}] + 
       Tensor[rm,{i,a[j],k,b[j]},{l,x}] /; RiemannQ[rm],
  HoldPattern[Tensor[rm_,{k_,l_,a_[i_],j_},{b_[i_],x___}]] :>
       Tensor[rm,{a[i],j,b[i],l},{k,x}] + 
       Tensor[rm,{a[i],j,k,b[i]},{l,x}] /; RiemannQ[rm],
  HoldPattern[Tensor[rm_,{k_,l_,i_,a_[j_]},{b_[j_],x___}]] :>
       Tensor[rm,{i,a[j],b[j],l},{k,x}] + 
       Tensor[rm,{i,a[j],k,b[j]},{l,x}] /; RiemannQ[rm],
  HoldPattern[Tensor[rc_,{a_[i_],j_},{b_[i_],x___}]] :>
       (1/2) Tensor[rc,{a[i],b[i]},{j,x}] /; RiemannQ[rc],
  HoldPattern[Tensor[rc_,{i_,a_[j_]},{b_[j_],x___}]] :>
       (1/2) Tensor[rc,{a[j],b[j]},{i,x}] /; RiemannQ[rc]
  };

(* BianchiRules allows the user to selectively apply the first or second
   Bianchi identity on specific indices.  The first Bianchi identity is
   applied to any tensor with Symmetries -> RiemannSymmetries; the second
   only to tensors that have been defined as Riemannian curvature tensors. *)

BianchiRules[i_,j_,k_] := {
  HoldPattern[Tensor[rm_,{i,j,k,l_},{x___}]] :> 
    -Tensor[rm,{j,k,i,l},{x}] - Tensor[rm,{k,i,j,l},{x}] /; 
    (Symmetries /. TensorData[rm]) === RiemannSymmetries,
  HoldPattern[Tensor[rm_,{i,j,l_,k},{x___}]] :> 
    Tensor[rm,{j,k,i,l},{x}] + Tensor[rm,{k,i,j,l},{x}] /; 
    (Symmetries /. TensorData[rm]) === RiemannSymmetries,
  HoldPattern[Tensor[rm_,{m_,i,j,k},{x___}]] :> 
    -Tensor[rm,{m,j,k,i},{x}] - Tensor[rm,{m,k,i,j},{x}] /; 
    (Symmetries /. TensorData[rm]) === RiemannSymmetries,
  HoldPattern[Tensor[rm_,{i,m_,j,k},{x___}]] :> 
    Tensor[rm,{m,j,k,i},{x}] + Tensor[rm,{m,k,i,j},{x}] /; 
    (Symmetries /. TensorData[rm]) === RiemannSymmetries,
  HoldPattern[Tensor[rm_,{i,j,p_,l_},{k,x___}]] :>
    -Tensor[rm,{k,i,p,l},{j,x}] - Tensor[rm,{j,k,p,l},{i,x}] /; RiemannQ[rm],
  HoldPattern[Tensor[rm_,{p_,l_,i,j},{k,x___}]] :>
    -Tensor[rm,{p,l,k,i},{j,x}] - Tensor[rm,{p,l,j,k},{i,x}] /; RiemannQ[rm]
  };

(********************** CurvToConnRule ********************************)

CurvToConnRule := {
  HoldPattern[ (r:Tensor[ rm_, {_,_,_,_}, {___} ]) /; RiemannQ[rm] ] :>
    (LowerAllIndices[r] /.
     HoldPattern[Tensor[rm,{L[i_],L[j_],L[k_],L[l_]},{x___}]] :>
      If[ RiemannConvention[Bundle[i]] === SecondUp,
          (*then*) 1,
          (*else*) -1] * 
      expandedCovD[ 
          (Conn[L[i], U[#], L[l]]*Conn[L[#], L[j], L[k]] - 
           Conn[L[i], U[#], L[k]]*Conn[L[#], L[j], L[l]] + 
           Del[Basis[L[k]], Conn[L[i], U[#], L[l]]] *
             Metric[Bundle[j]] [L[j], L[#]] - 
           Del[Basis[L[l]], Conn[L[i], U[#], L[k]]] *
             Metric[Bundle[j]] [L[j], L[#]] &) @
           NewBundleSymbol[Bundle[i]] +
          (Conn[L[i], L[j], L[#]]*Conn[L[k], U[#], L[l]] &) @
           NewBundleSymbol[Bundle[k]] - 
          (Conn[L[i], L[j], L[#]]*Conn[L[l], U[#], L[k]] &) @
           NewBundleSymbol[Bundle[l]],
        {x} ]),
  HoldPattern[ (r:Tensor[ rc_, {_,_}, {___} ]) /; RiemannQ[rc] ] :>
    (LowerAllIndices[r] /.
     HoldPattern[Tensor[rc,{L[i_],L[l_]},{x___}]] :>
      expandedCovD[ 
          Del[Basis[L[#1]], Conn[L[i], U[#1], L[l]]] - 
          Del[Basis[L[l]], Conn[L[i], U[#1], L[#1]]] + 
          Conn[L[i], U[#1], L[l]]*Conn[L[#1], U[#2], L[#2]] - 
          Conn[L[i], U[#1], L[#2]]*Conn[L[l], U[#2], L[#1]] & @@
          Table[NewBundleSymbol[Bundle[i]],{2}],
        {x} ]),
  HoldPattern[ (s:Tensor[ sc_, {}, {___} ]) /; RiemannQ[sc] ] :>
    (LowerAllIndices[s] /.
     HoldPattern[Tensor[sc,{},{x___}]] :> Module[{g},
      g = Riemannian /. TensorData[sc];
      expandedCovD[ 
          Conn[L[#1], U[#2], L[#2]]*Conn[L[#3], U[#1], U[#3]] - 
          Conn[L[#1], U[#2], L[#3]]*Conn[U[#1], U[#3], L[#2]] - 
          Del[Basis[L[#1]], Conn[L[#3], U[#2], L[#2]]]*g[U[#1], U[#3]] + 
          Del[Basis[L[#1]], Conn[L[#2], U[#1], L[#3]]]*g[U[#2], U[#3]] & @@
          Table[ NewBundleSymbol[Bundles[g][[1,1]]], {3}],
        {x} ] ]),
  HoldPattern[ c:Tensor[ Curv, {_,_,_,_}, {___} ] ] :>
    (LowerAllIndices[c] /.
     HoldPattern[Tensor[Curv,{L[i_],L[j_],L[k_],L[l_]},{x___}]] :>
      expandedCovD[ 
          (Conn[L[i], U[#], L[l]]*Conn[L[#], L[j], L[k]] - 
           Conn[L[i], U[#], L[k]]*Conn[L[#], L[j], L[l]] + 
           Del[Basis[L[k]], Conn[L[i], U[#], L[l]]] *
             Metric[Bundle[j]] [L[j], L[#]] - 
           Del[Basis[L[l]], Conn[L[i], U[#], L[k]]] *
             Metric[Bundle[j]] [L[j], L[#]] &) @
           NewBundleSymbol[Bundle[i]] +
          (Conn[L[i], L[j], L[#]]*Conn[L[k], U[#], L[l]] &) @
           NewBundleSymbol[Bundle[k]] - 
          (Conn[L[i], L[j], L[#]]*Conn[L[l], U[#], L[k]] &) @
           NewBundleSymbol[Bundle[l]],
        {x} ])
  };

(*********************** ConnToMetricRule **********************************)

ConnToMetricRule := {
  HoldPattern[Tensor[Conn,{L[i_],j_,L[k_]},{}]] :> Module[ {g},
    g = Metric[Bundle[j]];
      (1/2) g[j,U[#]] * (
          Del[ Basis[L[i]], g[L[k],L[#]] ] +
          Del[ Basis[L[k]], g[L[i],L[#]] ] -
          Del[ Basis[L[#]], g[L[i],L[k]] ]) & @
      NewBundleSymbol[Conjugate[Bundle[j]]]
    ] /; MatchQ[ Union[{Bundle[L[i]],Bundle[j],Bundle[L[k]]}],
                 {b_} | {b_,Conjugate[b_]} | {Conjugate[b_],b_} ] &&
         CommutingFrame[Bundle[j]],
  HoldPattern[Tensor[Conn,{U[i_],j_,k_},{}]] :> 
    (Tensor[Conn,{U[i],j,k},{}] /. LowerIndicesRule[U[i]] /.
                                   ConnToMetricRule) /; 
    MatchQ[ Union[{Bundle[U[i]],Bundle[j],Bundle[k]}],
                 {b_} | {b_,Conjugate[b_]} | {Conjugate[b_],b_} ] &&
         CommutingFrame[Bundle[j]], 
  HoldPattern[Tensor[Conn,{i_,j_,U[k_]},{}]] :> 
    (Tensor[Conn,{i,j,U[k]},{}] /. LowerIndicesRule[U[k]] /.
                                   ConnToMetricRule) /;
    MatchQ[ Union[{Bundle[i],Bundle[j],Bundle[U[k]]}],
                 {b_} | {b_,Conjugate[b_]} | {Conjugate[b_],b_} ] &&
         CommutingFrame[Bundle[j]]
  };    

(*********************** Curvatures for arbitrary metrics ******************)

RiemannTensor[h_] := (
  Print["ERROR: ",h," is not a symmetric 2-tensor"];
  ERROR[RiemannTensor[h]]) /; !Rank[h]===2 || !SymmetricQ[h];

HoldPattern[RiemannTensor[ f_. (g_?MetricQ) ] /; ConstantQ[f] ] :=
  f RiemannTensor[Bundles[g][[1,1]]];

RiemannTensor/: Conjugate[HoldPattern[RiemannTensor[h_]]] := RiemannTensor[h];

RiemannTensor/: TCompute[ HoldPattern[RiemannTensor[h_]], 
                          {i_,j_,k_,l_} ] := Module[{b,rmsign,dums},
  b = Inverse[h];
  rmsign = If[RiemannConvention[ Bundles[h][[1,1]] ] === SecondUp,
              (*then*)   1,
              (*else*)   -1];
  dums = Table[NewBundleSymbol /@ (Bundles[h] [[1]] ), {2}];
  Return[ rmsign * (
    (1/2) h[i, k] [l, j] - 
    (1/2) h[i, l] [k, j] - 
    (1/2) h[j, k] [l, i] + 
    (1/2) h[j, l] [k, i] + 
    Plus @@ Flatten[ Outer @@ Join[ {
        (1/2) Curv[i, j, k, L[#1]] h[l, U[#1]] - 
        (1/2) Curv[i, j, l, L[#1]] h[k, U[#1]] &},
      Take[dums,1] ] ] + 
    Plus @@ Flatten[ Outer @@ Join[ {
        (1/4) h[i, k] [L[#1]] h[j, l] [L[#2]] b[U[#1], U[#2]] - 
        (1/4) h[i, L[#1]] [k] h[j, l] [L[#2]] b[U[#1], U[#2]] - 
        (1/4) h[i, l] [L[#1]] h[j, k] [L[#2]] b[U[#1], U[#2]] + 
        (1/4) h[i, L[#1]] [l] h[j, k] [L[#2]] b[U[#1], U[#2]] + 
        (1/4) h[i, l] [L[#1]] h[j, L[#2]] [k] b[U[#1], U[#2]] - 
        (1/4) h[i, L[#1]] [l] h[j, L[#2]] [k] b[U[#1], U[#2]] - 
        (1/4) h[i, k] [L[#1]] h[j, L[#2]] [l] b[U[#1], U[#2]] + 
        (1/4) h[i, L[#1]] [k] h[j, L[#2]] [l] b[U[#1], U[#2]] - 
        (1/4) h[j, l] [L[#1]] h[k, L[#2]] [i] b[U[#1], U[#2]] + 
        (1/4) h[j, L[#1]] [l] h[k, L[#2]] [i] b[U[#1], U[#2]] + 
        (1/4) h[i, l] [L[#1]] h[k, L[#2]] [j] b[U[#1], U[#2]] - 
        (1/4) h[i, L[#1]] [l] h[k, L[#2]] [j] b[U[#1], U[#2]] + 
        (1/4) h[j, k] [L[#1]] h[l, L[#2]] [i] b[U[#1], U[#2]] - 
        (1/4) h[j, L[#1]] [k] h[l, L[#2]] [i] b[U[#1], U[#2]] - 
        (1/4) h[k, L[#1]] [j] h[l, L[#2]] [i] b[U[#1], U[#2]] - 
        (1/4) h[i, k] [L[#1]] h[l, L[#2]] [j] b[U[#1], U[#2]] + 
        (1/4) h[i, L[#1]] [k] h[l, L[#2]] [j] b[U[#1], U[#2]] + 
        (1/4) h[k, L[#1]] [i] h[l, L[#2]] [j] b[U[#1], U[#2]] &},
      dums ] ] 
   )];
  ];

RicciTensor[h_] := (
  Print["ERROR: ",h," is not a symmetric 2-tensor"];
  ERROR[RicciTensor[h]]) /; !Rank[h]===2 || !SymmetricQ[h];

HoldPattern[RicciTensor[f_. * (g_?MetricQ)] /; ConstantQ[f]] := 
  RicciTensor[Bundles[g][[1,1]]];

RicciTensor/: Conjugate[HoldPattern[RicciTensor[h_]]] := RicciTensor[h];

RicciTensor/: TCompute[ HoldPattern[RicciTensor[h_]], {i_,j_} ] := 
  Module[{b,dums},
  b = Inverse[h];
  dums = Table[NewBundleSymbol /@ (Bundles[h] [[1]] ), {4}];
  Return[ 
    Plus @@ Flatten[ Outer @@ Join[ {
        (1/2) Curv[i, L[#1], U[#1], j] &},
      Take[dums,1] ] ] + 
    Plus @@ Flatten[ Outer @@ Join[ {
        (1/2) Curv[i, L[#1], L[#2], L[#3]] h[j, U[#3]] b[U[#1], U[#2]] &},
      Take[dums,3] ] ] +
    Plus @@ Flatten[ Outer @@ Join[ { 
        (1/2) h[i, L[#1]] [j, L[#2]] b[U[#1], U[#2]] - 
        (1/2) h[i, j] [L[#1], L[#2]] b[U[#1], U[#2]] - 
        (1/2) h[L[#1], L[#2]] [j, i] b[U[#1], U[#2]] + 
        (1/2) h[L[#1], j] [L[#2], i] b[U[#1], U[#2]] &},
      Take[dums,2] ] ] +
    Plus @@ Flatten[ Outer @@ Join[ { 
        (1/4) h[i, L[#1]] [j] h[L[#2], L[#3]] [L[#4]] *
              b[U[#1], U[#4]] * b[U[#2], U[#3]] - 
        (1/2) h[i, L[#1]] [j] h[L[#2], L[#3]] [L[#4]] *
              b[U[#1], U[#2]] * b[U[#3], U[#4]] - 
        (1/4) h[i, j] [L[#1]] h[L[#2], L[#3]] [L[#4]] *
              b[U[#1], U[#4]] * b[U[#2], U[#3]] + 
        (1/2) h[i, j] [L[#1]] h[L[#2], L[#3]] [L[#4]] *
              b[U[#1], U[#2]] * b[U[#3], U[#4]] - 
        (1/2) h[L[#1], j] [i] h[L[#2], L[#3]] [L[#4]] *
              b[U[#1], U[#2]] * b[U[#3], U[#4]] - 
        (1/2) h[i, L[#1]] [L[#2]] h[j, L[#3]] [L[#4]] *
              b[U[#1], U[#4]] * b[U[#2], U[#3]] + 
        (1/2) h[i, L[#1]] [L[#2]] h[j, L[#3]] [L[#4]] *
              b[U[#1], U[#3]] * b[U[#2], U[#4]] - 
        (1/4) h[L[#1], L[#2]] [i] h[j, L[#3]] [L[#4]] *
              b[U[#1], U[#4]] * b[U[#2], U[#3]] + 
        (1/4) h[L[#1], L[#2]] [L[#3]] h[j, L[#4]] [i] *
              b[U[#1], U[#2]] * b[U[#3], U[#4]] + 
        (1/4) h[L[#1], j] [L[#2]] h[L[#3], L[#4]] [i] * 
              b[U[#1], U[#3]] * b[U[#2], U[#4]] + 
        (1/4) h[L[#1], L[#2]] [i] h[L[#3], L[#4]] [j]  *
              b[U[#1], U[#3]] * b[U[#2], U[#4]]
        &},
      dums
  ]]];
  ];

ScalarCurv[h_] := (
  Print["ERROR: ",h," is not a symmetric 2-tensor"];
  ERROR[ScalarCurv[h]]) /; !Rank[h]===2 || !SymmetricQ[h];

HoldPattern[ScalarCurv[f_. * (g_?MetricQ)] /; ConstantQ[f]] := 
  f^(-1) ScalarCurv[Bundles[g][[1,1]]];

ScalarCurv/: Conjugate[HoldPattern[ScalarCurv[h_]]] := ScalarCurv[h];

ScalarCurv/: TCompute[ HoldPattern[ScalarCurv[h_]], {} ] := 
  Module[{b,ricci,dums},
  b = Inverse[h];
  dums = Table[NewBundleSymbol /@ (Bundles[h] [[1]] ), {6}];
  Return[ 
  Plus @@ Flatten[ Outer @@ Join[ {
      Curv[L[#1], L[#2], U[#2], L[#3]] b[U[#1], U[#3]] &},
    Take[dums,3] ] ] + 
  Plus @@ Flatten[ Outer @@ Join[ {
      h[L[#1], L[#2]] [L[#3], L[#4]] b[U[#1], U[#3]] * b[U[#2], U[#4]] - 
      h[L[#1], L[#2]] [L[#3], L[#4]] b[U[#1], U[#2]] * b[U[#3], U[#4]] &},
    Take[dums,4] ] ] +
  Plus @@ Flatten[ Outer @@ Join[ { 
    - (1/2) h[L[#1], L[#2]] [L[#3]] h[L[#4], L[#5]] [L[#6]] *
            b[U[#1], U[#4]] * b[U[#2], U[#6]] * b[U[#3], U[#5]] + 
      (3/4) h[L[#1], L[#2]] [L[#3]] h[L[#4], L[#5]] [L[#6]] *
            b[U[#1], U[#4]] * b[U[#2], U[#5]] * b[U[#3], U[#6]] + 
      (1/2) h[L[#1], L[#2]] [L[#3]] h[L[#4], L[#5]] [L[#6]] *
            b[U[#1], U[#3]] * b[U[#2], U[#6]] * b[U[#4], U[#5]] - 
      (1/4) h[L[#1], L[#2]] [L[#3]] h[L[#4], L[#5]] [L[#6]] *
            b[U[#1], U[#2]] * b[U[#3], U[#6]] * b[U[#4], U[#5]] - 
      h[L[#1], L[#2]] [L[#3]] h[L[#4], L[#5]] [L[#6]] *
            b[U[#1], U[#3]] * b[U[#2], U[#4]] * b[U[#5], U[#6]] + 
      (1/2) h[L[#1], L[#2]] [L[#3]] h[L[#4], L[#5]] [L[#6]] *
            b[U[#1], U[#2]] * b[U[#3], U[#4]] * b[U[#5], U[#6]]
      &},
    dums
  ]]];
  ];

LeviCivitaConnection[h_] := (
  Print["ERROR: ",h," is not a symmetric 2-tensor"];
  ERROR[LeviCivitaConnection[h]]) /; !Rank[h]===2 || !SymmetricQ[h];

LeviCivitaConnection/: TCompute[ HoldPattern[LeviCivitaConnection[h_]], 
                                 {i_,j_,k_} ] := 
  InsertIndices[ Conn, {i,j,k} ] +
  Plus @@ Flatten[ Outer @@ Join[ {
      (1/2) InsertIndices[ Inverse[h], {j,U[#]} ] *
              CovD[ InsertIndices[ h, {i,L[#]} ], {k} ] +
      (1/2) InsertIndices[ Inverse[h], {j,U[#]} ] *
              CovD[ InsertIndices[ h, {k,L[#]} ], {i} ] -
      (1/2) InsertIndices[ Inverse[h], {j,U[#]} ] *
              CovD[ InsertIndices[ h, {i,k} ], {L[#]} ]
    &},
    {NewBundleSymbol /@ (Bundles[h] [[1]] )}
  ]];

LeviCivitaConnection/: TCompute[ HoldPattern[LeviCivitaConnection[h_]],
                                 {i_,j_} ] := 
  InsertIndices[ Conn, {i,j} ] +
  Plus @@ Flatten[ Outer @@ Join[ {
      (1/2) InsertIndices[ Inverse[h], {j,U[#1]} ] *
             CovD[ InsertIndices[ h, {i,L[#1]} ], {L[#2]} ] * Basis[U[#2]] +
      (1/2) InsertIndices[ Inverse[h], {j,U[#1]} ] *
             CovD[ InsertIndices[ h, {L[#2],L[#1]} ], {i} ] * Basis[U[#2]] -
      (1/2) InsertIndices[ Inverse[h], {j,U[#1]} ] *
             CovD[ InsertIndices[ h, {i,L[#2]} ], {L[#1]} ] * Basis[U[#2]]
    &},
    Table[ NewBundleSymbol /@ (Bundles[h] [[1]] ), {2} ]
  ]];

LeviCivitaConnection/: TCompute[ HoldPattern[LeviCivitaConnection[h_]],
                                 {i___} ] := (
  Message[ Index::error, "LeviCivitaConnection" ];
  Return[ERROR[ InsertIndices[ LeviCivitaConnection[h], {i} ]]];
  ) /; Length[{i}] =!= 2 && Length[{i}] =!= 3;


(*************************** Tensor *****************************)

SetAttributes[{slotsfullQ,decompOK,bundleOK,nonscalarQ}, HoldAll];

(***************** Tensor index-eating behavior **************)

Tensor/: TCompute[ HoldPattern[Tensor[t_,{i___},{j__}]], {k___} ] :=
  Tensor[t,{i},{j,k}];

Tensor/: TCompute[ HoldPattern[Tensor[t_,{i___},{}]], {j___}] := Which[
  slotsfullQ[Tensor[t,{i},{}]],
      Tensor[t,{i},{j}],
  !decompOK[Tensor[t,{i},{}],{j}],
      ERROR[Tensor[t,{i},{}][j]],
  True,
      Tensor[t,{i,j},{}]
  ];

(* decompOK[ tensorname, {existing indices}, {new indices} ] returns True
   if the number of new indices provided is consistent with the rank
   of tensorname. *)

decompOK[Tensor[Inverse[_],{i___},{}],{j___}] := {i}==={} && Length[{j}]===2;

decompOK[Tensor[t_,{i___},{}],{j___}] := Module[{compressedrank},
  Switch[Head[TensorRankList[t]],
    Integer,
      If[{i}==={} && Length[{j}]===TensorRankList[t],
         (*then*) Return[True],
         (*else*) Message[Index::error,Tensor[t,{i},{}]];
                  Return[False]
        ],
    List,
      (* restructure list of ranks by consolidating ranks of already filled slots. *)
      compressedrank = TensorRankList[t] //.
        {{a_,b___} :> {0,a,b} /; Length[{i}] == 0 && a > 0,
         {a_,b_,c___} :> {a+b,c} /; a+b <= Length[{i}],
         {a_,b_,c_,d___} :> {a,b+c,d} /; a >= Length[{i}] && b+c <= Length[{j}]
        };
      If[Length[compressedrank] >= 2 &&
          compressedrank[[2]]===Length[{j}] && compressedrank[[1]]===Length[{i}],
        (*then*) Return[True],
        (*else*) Message[Index::error,Tensor[t,{i},{}]];
                 Return[False]
        ]
    ];
  ];

(**************** Conjugation rules for tensors **********************)

Tensor/: Conjugate[Tensor[Inverse[x_],i_,j_]] :=
         Tensor[Inverse[Conjugate[x]],Conjugate[i],Conjugate[j]];

Tensor/: Conjugate[Tensor[t_[Bar],i_,j_]] :=
         Tensor[t,Conjugate[i],Conjugate[j]];

Tensor/: Conjugate[Tensor[t_,i_,j_]] :=
         Switch[Type[t],
                {___,Real,___},  Tensor[t,Conjugate[i],Conjugate[j]],
                {Imaginary},     -Tensor[t,Conjugate[i],Conjugate[j]],
                {Complex},       Tensor[t[Bar],Conjugate[i],Conjugate[j]],
                _,               ERROR[Conjugate[Tensor[t,i,j]]]
         ];

(******************** Inverse tensors *********************************)

(* Inverse has to be handled a little differently from the other Ricci
   functions when indices are inserted, because there's no way to
   express the components of the inverse in terms of the components of
   the original tensor using our index conventions.  In effect, we
   manufacture a new tensor whose name is "Inverse[x]". *)

(* An indexed tensor times its inverse gives the identity *)

Tensor/: Tensor[Inverse[Tensor[x_,{},{}]], {i_,a_[j_]}, {}] *
         Tensor[x_, {b_[j_],k_}, {}] :=
           Tensor[Kronecker,{i,k},{}] /; PairQ[a[j],b[j]];

Tensor/: Tensor[Inverse[Tensor[x_,{},{}]], {a_[j_],i_}, {}] *
         Tensor[x_, {k_,b_[j_]}, {}] :=
           Tensor[Kronecker,{k,i},{}] /; PairQ[a[j],b[j]];

Tensor/: Tensor[Inverse[Tensor[x_,{},{}]], {i___,a_[j_],k___}, {}] *
         Tensor[x_, {l___,b_[j_],m___}, {}] :=
           Tensor[Kronecker,{i,k,l,m},{}] /; SymmetricQ[x] && PairQ[a[j],b[j]];

Tensor/: Tensor[Inverse[Tensor[x_,{},{}]], {i___,a_[j_],k___}, {}] *
         Tensor[x_, {l___,b_[j_],m___}, {}] :=
           - Tensor[Kronecker,{i,k,l,m},{}] /; SkewQ[x] && PairQ[a[j],b[j]];

(*********** Functions for determining tensor properties **************)

TensorData[t_[Bar]] := TensorData[t];

TensorRankList[Tensor[t_,_,_]] := Rank /. TensorData[t];
TensorRankList[t_[Bar]] := TensorRankList[t];

Type[Tensor[Inverse[x_],_,_]] := {Real} /; x===Conjugate[x];
Type[Tensor[t_[Bar],i_,j_]] := Type[Tensor[t,{},{}]];
Type[Tensor[t_,_,_]] := Type /. TensorData[t];

RiemannQ[Tensor[t_[Bar],_,_]] := False;
RiemannQ[Tensor[Inverse[x_],_,_]] := False;
RiemannQ[Tensor[t_,_,_]] :=
  (Riemannian /. TensorData[t]) =!= False;

Unprotect[TensorQ];
TensorQ[_] := False;
TensorQ[t_[Bar]] := TensorQ[t];
Protect[TensorQ];

slotsfullQ[Tensor[Inverse[_],{i___},{___}]] :=
  Length[Unevaluated[i]]===2;
slotsfullQ[Tensor[t_[Bar],{i___},{___}]] :=
  TotalRank[t]===Length[Unevaluated[i]];
slotsfullQ[Tensor[t_,{i___},{___}]] :=
  TotalRank[t]===Length[Unevaluated[i]];

nonscalarQ[t_] := !slotsfullQ[t];

(********** Argument-checking subroutines for DefineTensor ***********)

symQ[{s__},{d__}] := And @@ (symQ @@ # &) /@ Transpose[{{s},{d}}]  /;
  Length[{s}] === Length[{d}] &&
  !MemberQ[ {s}, RiemannSymmetries | _TensorSymmetry ];
symQ[s_,{d__}] := SameQ[s, NoSymmetries];
symQ[RiemannSymmetries,d_] := d===4;
symQ[TensorSymmetry[_,d_,_],d_] := TrueQ[d > 1];
symQ[Hermitian,d_] := d===2;
symQ[SkewHermitian,d_] := d===2;
symQ[s_,_] := (s===Symmetric || s===Alternating ||
               s===NoSymmetries);
symQ[___] = False;

Skew = Alternating;

tentypeQ[t_List,d_] :=
  (t==={Complex}) ||
  (t==={Imaginary}) ||
  (t==={Real}) ||
  (d===0 &&
   ContainsQ[ {Real,Positive,Negative,NonPositive,NonNegative}, t ] &&
   !ContainsQ[ t, {Positive,Negative} ] &&
   !ContainsQ[ t, {Positive,NonPositive} ] &&
   !ContainsQ[ t, {Negative,NonNegative} ]);
tentypeQ[t_,d_] := tentypeQ[{t},d];

degQ[{d__}] := And @@ (degQ /@ {d});
degQ[d_] := IntegerQ[d] && d >= 0;

bundleQ[{s__},{d__}] := And @@ (bundleQ @@ # &) /@ Transpose[{{s},{d}}]  /;
  Length[{s}] === Length[{d}];
bundleQ[b_,{__}] := bundleQ[b];
bundleQ[{b__},d_Integer] := And @@ bundleQ /@ {b};
bundleQ[b_,_] := bundleQ[b];
bundleQ[Automatic] := True;
bundleQ[Same] := True;
bundleQ[Any] := True;
bundleQ[Null] := True;
bundleQ[b_] := BundleQ[b];
bundleQ[___] := False;

varQ[Covariant,_] := True;
varQ[Contravariant,_] := True;
varQ[{v___},d_] := Length[{v}]===(Plus @@ Flatten[{d}]) &&
                   And @@ (varQ[#,0]&) /@ {v};
varQ[Same,{_,__}] := True;
varQ[___] := False;

Co = Covariant;
Con = Contravariant;

Any/:  Conjugate[Any]  = Any;
Same/: Conjugate[Same] = Same;

(************************** DefineTensor **************************)

Options[DefineTensor] := {Symmetries -> NoSymmetries,
                         Type -> Real,
                         TeXFormat -> "",
                         Quiet -> $Quiet,
                         Bundle -> Null,
                         Variance -> Covariant,
                         Riemannian -> False};

DefineTensor[] := Message[DefineTensor::argm,"DefineTensor",0,2];
DefineTensor[_] := Message[DefineTensor::argmu,"DefineTensor",2];

DefineTensor[l_List, d_, opts___] := DefineTensor[#,d,opts]& /@ l;

DefineTensor[t_, deg_, opts___] := Module[{sym,type,tex,quiet,bun,var,riem},
   If [!optionsOK[DefineTensor,opts], Return[]];
   {sym,type,tex,quiet,bun,var,riem} =
     {Symmetries,Type,TeXFormat,Quiet,Bundle,Variance,Riemannian} /.
     {opts} /. Options[DefineTensor];
   If[ bun===Null,
     bun =
       If[Head[deg]===List,
         (*then*) Table[$DefaultTangentBundle, {Length[deg]}],
         (*else*) $DefaultTangentBundle]
     ];
   Which[
    !newsymQ[t],
      Message[DefineTensor::invalid,HoldForm[t],"tensor name"],
    !degQ[deg],
      Message[DefineTensor::invalid,deg,"tensor degree"],
    !symQ[sym,deg],
      Message[DefineTensor::invalid,sym,"tensor symmetry for "<>
                                         ToString[HoldForm[t]]],
    !tentypeQ[type,deg],
      Message[DefineTensor::invalid,type,"tensor type for "<>
                                         ToString[HoldForm[t]]],
    !SameQ[Head[tex],String],
      Message[DefineTensor::invalid,tex,"tensor TeX format"],
    !FreeQ[bun,Null],
      Message[Bundle::error,t],
    !bundleQ[bun,deg],
      Message[DefineTensor::invalid,bun,
              "bundle name for "<>ToString[HoldForm[t]]],
    !varQ[var,deg],
      Message[DefineTensor::invalid,var,
              "variance for "<>ToString[HoldForm[t]]],
    !tfQ[riem],
      Message[DefineTensor::invalid,riem,"value for the Riemannian option"],
    (* the arguments pass, so define the thing *)
    True,
      (* for (skew-) Hermitian tensors, force Type = Real *)
      If[sym==Hermitian || sym==SkewHermitian, type = Real];
      t = Tensor[t,{},{}];
      t/: TensorQ[t] = True;
      t/: TotalRank[t] = Plus @@ Flatten[{deg}];
      t/: TensorData[t] = {Rank       -> deg,
                           Symmetries -> sym,
                           Type       -> Null,  (*Declare will change later*)
                           Bundles    -> Null,  (*ditto*)
                           TangentBundle -> TangentBundle[
                                            Flatten[{bun,Null}] [[1]] /.
                                            {Automatic->Null,
                                            Same->Null,
                                            Any->Null}],
                           Variance   -> If[Head[var]===List,
                                         (*then*) var,
                                         (*else*) Table[var,{TotalRank[t]}]],
                           TeXFormat  -> tex,
                           Riemannian -> If[ riem,
                                             Metric[Flatten[{bun}][[1]]],
                                             False ]};
      defSymmetryRelations[t,sym];
      Protect[t];
      DeclareTensor[t,Type->type,Bundle->bun,TeXFormat->tex];
      If[!quiet,Print["Tensor ", HoldForm[t], " defined."];
                Print["  Rank = ", TensorRankList[t],
                      "  Symmetries = ", Symmetries /. TensorData[t]];
                Print["  Type = ", Type[t], "  Bundle = ",bun,
                      "  Variance = ", var]
        ];
      Return[t]
  ]];

(* The following function turns the user's Bundle option into a list,
   each entry of which is a sublist representing the allowed bundles
   for that index slot. *)

makeBundleList[Tensor[t_,{},{}],buns_,ranklist_] := Module[{bunlist},
  bunlist =
  If[Head[buns]===List && Head[ranklist]===List,
    (*then there are lists of bundles & degrees*)
      Flatten[ (Table @@ # &) /@
          Transpose[ {Flatten /@ List /@ buns, List /@ ranklist} ], 1
          ],
    (*else there's only one bundle, or one list to be repeated*)
      Table[ Flatten[{buns}], {TotalRank[t]} ]
    ];
  (* If the tensor is real or imaginary, allow conjugates of specified bundles too. *)
  If[ Type[t]=!={Complex},
    (*then*) bunlist = (Union[Flatten[{#, Conjugate/@#}]]&) /@ bunlist];
  Return[bunlist];
  ];

(********************** DeclareTensor *************************)

(* Called by Declare when the user gives Declare a tensor name
   as first argument.  (See Constant.m) *)

Options[DeclareTensor] := {Type      -> Null,
                           TeXFormat -> Null,
                           Bundle    -> Null};

DeclareTensor[Tensor[t_,{},{}], opts___ ] := Module[{type,bun,tex},
  {type,bun,tex} =
    ({Type,Bundle,TeXFormat} /. {opts} /. Options[DeclareTensor]);
  Which[
    type =!= Null && !tentypeQ[type,TensorRankList[t]],
      Message[Declare::invalid, type, "tensor type for "<>
                                      ToString[HoldForm[t]]],
    bun =!= Null && !bundleQ[bun,TensorRankList[t]],
      Message[Declare::invalid, bun,
              "bundle name for "<>ToString[HoldForm[t]]],
    True,
      Unprotect[t];
      If[ type =!= Null,
          type = Flatten[{type}];
          If[ MemberQ[type, Odd|Even|Integer|Positive|
                            Negative|NonPositive|
                            NonNegative],
              (*then*) type = Union[type,{Real}] ];
          If[ MemberQ[type, Odd|Even],
              (*then*) type = Union[type,{Integer}] ];
          t/: TensorData[t] = TensorData[t] /.
                          {HoldPattern[Type -> _] :> (Type -> type)};
        ];
      If[ bun =!= Null,
          t/: TensorData[t] = TensorData[t] /.
                          {HoldPattern[Bundles -> _] :>
                           (Bundles -> 
                             makeBundleList[t,bun,TensorRankList[t]])}
        ];
      If[ tex =!= Null,
          t/: TensorData[t] = TensorData[t] /.
                          {HoldPattern[TeXFormat -> _] :>
                           (TeXFormat -> tex)}
        ];
      If[ Type[t]=!={Complex},
          t/:  TensorData[t] = TensorData[t] /.
                   {HoldPattern[Bundles -> b_] :>
                    (Bundles -> ((Union[Flatten[{#,
                                  Conjugate/@#}]]&) /@ b))}
        ];
    ];
  Return[t];
  ];

(********************** UndefineTensor *************************)

Attributes[UndefineTensor] = {Listable};

Options[UndefineTensor] := {Quiet -> $Quiet};

UndefineTensor[Tensor[t_,{},{}],opts___] := Module[{quiet},
  If[ !TensorQ[t],
        Message[UndefineTensor::invalid, HoldForm[t], "tensor name"];
        Return[]
    ];
  If[ !optionsOK[UndefineTensor,opts], Return[]];
  quiet = Quiet /. {opts} /. Options[UndefineTensor];
  If[ !quiet,Print["Undefining Tensor ",HoldForm[t]]];
  Unprotect[t];
  ClearAll[t];
  ];

UndefineTensor[t_,___] :=
  Message[UndefineTensor::invalid, HoldForm[t], "tensor name"];

(*********************** Symmetries ****************************)

(* If symmetries are specified as a list, then call "defSymmetryRelations"
   once for each entry in the list.  Each call sets up the symmetry
   relations for a subset of the index positions. *)

defSymmetryRelations[Tensor[t_,{},{}], x_List] := Module[{n,pos=0},
  Do[ (* for n=1 to Length[TensorRankList[t]] *)
    defSymmetryRelations[t, x[[n]], pos, TensorRankList[t][[n]] ];
    pos += TensorRankList[t][[n]],
    {n, Length[TensorRankList[t]]}
    ];
  ];

(* If symmetries = NoSymmetries, do nothing *)

defSymmetryRelations[_, NoSymmetries, ___] := Null;

(* If 3rd & 4th arguments p,q are present, they indicate that the symmetries
   are to apply only to the q index positions starting at position p+1. *)

defSymmetryRelations[Tensor[t_,{},{}], Symmetric] := (
  t/: HoldPattern[Tensor[t,{i_,j__},{x___}]] :=
         Tensor[t, IndexSort[{i,j}], {x}] /; !IndexOrderedQ[{i,j}];
  If[Type[t] === {Complex},
     t/: HoldPattern[Tensor[t[Bar],{i_,j__},{x___}]] :=
         Tensor[t[Bar], IndexSort[{i,j}], {x}] /; !IndexOrderedQ[{i,j}]]
  );

defSymmetryRelations[Tensor[t_,{},{}], Symmetric, p_, q_] := (
  t/: HoldPattern[Tensor[t,{i_,j__},{x___}]] :=
        Tensor[t,Join[Take[ {i,j}, p ],
                      IndexSort[ Take[ {i,j}, {p+1,p+q} ] ],
                      Drop[ {i,j}, p+q ]
                     ], {x}] /;
        Length[{i,j}] >= p+q &&
        !IndexOrderedQ[Take[{i,j},{p+1,p+q}]];
  If[Type[t]==={Complex},
     t/: HoldPattern[Tensor[t[Bar],{i_,j__},{x___}]] :=
         Tensor[t[Bar],Join[Take[ {i,j}, p ],
                            IndexSort[ Take[ {i,j}, {p+1,p+q} ] ],
                            Drop[ {i,j}, p+q ]
                           ], {x}] /;
         Length[{i,j}] >= p+q &&
         !IndexOrderedQ[Take[{i,j},{p+1,p+q}]]]
  );

defSymmetryRelations[Tensor[t_,{},{}], Alternating] := (
  t/: HoldPattern[Tensor[t,{i_,j__},{x___}]] :=
        IndexSignature[{i,j}] *
        Tensor[t,IndexSort[{i,j}],{x}] /; !IndexOrderedQ[{i,j}];
  t/: HoldPattern[Tensor[t,{___,L[a_],___,U[a_],___},{x___}]]:= 0 /;
        Type[Bundle[a]]===Real;
  t/: HoldPattern[Tensor[t,{___,a_,___,a_,___},{x___}]]:= 0;
  If[Type[t]==={Complex},
     t/: HoldPattern[Tensor[t[Bar],{i_,j__},{x___}]] :=
         IndexSignature[{i,j}] *
         Tensor[t[Bar],IndexSort[{i,j}],{x}] /; !IndexOrderedQ[{i,j}];
     t/: HoldPattern[Tensor[t[Bar],{___,L[a_],___,U[a_],___},{x___}]]:= 0 /;
         Type[Bundle[a]]===Real;
     t/: HoldPattern[Tensor[t[Bar],{___,a_,___,a_,___},{x___}]]:= 0 ]
  );

defSymmetryRelations[Tensor[t_,{},{}], Alternating, p_, q_] := (
  t/: HoldPattern[Tensor[t,{i_,j__},{x___}]] :=
        IndexSignature[Take[{i,j},{p+1,p+q}]] *
        Tensor[t,Join[Take[ {i,j}, p ],
                      IndexSort[ Take[{i,j}, {p+1,p+q}] ],
                      Drop[ {i,j}, p+q ]
                     ],
                 {x}] /;
        Length[{i,j}] >= p+q &&
        !IndexOrderedQ[Take[{i,j},{p+1,p+q}]];
  t/: HoldPattern[Tensor[t,{i___,L[a_],j___,U[a_],___},{___}]]:= 0 /;
        Type[Bundle[a]]===Real &&
        Length[{i}] >= p && Length[{i,L[a],j,U[a]}] <= p+q;
  t/: HoldPattern[Tensor[t,{i___,a_,j___,a_,___},{___}]]:= 0 /;
        Length[{i}] >= p && Length[{i,a,j,a}] <= p+q;
  If[Type[t]==={Complex},
     t/: HoldPattern[Tensor[t[Bar],{i_,j__},{x___}]] :=
           IndexSignature[Take[{i,j},{p+1,p+q}]] *
           Tensor[t[Bar],Join[Take[ {i,j}, p ],
                         IndexSort[ Take[{i,j}, {p+1,p+q}] ],
                         Drop[ {i,j}, p+q ]
                        ],
                    {x}] /;
           Length[{i,j}] >= p+q &&
           !IndexOrderedQ[Take[{i,j},{p+1,p+q}]];
     t/: HoldPattern[Tensor[t[Bar],{i___,L[a_],j___,U[a_],___},{___}]]:= 0 /;
           Type[Bundle[a]]===Real &&
           Length[{i}] >= p && Length[{i,L[a],j,U[a]}] <= p+q;
     t/: HoldPattern[Tensor[t[Bar],{i___,a_,j___,a_,___},{___}]]:= 0 /;
           Length[{i}] >= p && Length[{i,a,j,a}] <= p+q ]
    );

defSymmetryRelations[Tensor[t_,{},{}], Hermitian] := (
  defSymmetryRelations[ t, Symmetric ];
  MakeHermitian[ t ];
  );

defSymmetryRelations[Tensor[t_,{},{}], Hermitian, p_, q_] := (
  defSymmetryRelations[ t, Symmetric, p, q ];
  MakeHermitian[ t, p, q ];
  );

defSymmetryRelations[Tensor[t_,{},{}], SkewHermitian] := (
  defSymmetryRelations[ t, Alternating ];
  MakeHermitian[ t ];
  );

defSymmetryRelations[Tensor[t_,{},{}], SkewHermitian, p_, q_] := (
  defSymmetryRelations[ t, Alternating, p, q ];
  MakeHermitian[ t, p, q ];
  );

MakeHermitian[Tensor[t_,{},{}]] := (
  t/: HoldPattern[Tensor[t,{i_,j_},{___}]] := 0 /;
       SameQ @@ Bundle /@ Lower /@ {i,j};
  );

MakeHermitian[Tensor[t_,{},{}], p_, q_] := (
  t/: HoldPattern[Tensor[t,{i_,j__},{___}]] := 0 /;
        Length[{i,j}] >= p+q &&
        SameQ @@ Bundle /@ Lower /@ Take[ {i,j}, {p+1,p+q} ];
  );

defSymmetryRelations[Tensor[t_,{},{}],RiemannSymmetries] := (

  (** RULE 1
      This rule uses the symmetry

      T        = -T
       i j k l     j i k l

      to order the first two indices.  **)

  t/: HoldPattern[Tensor[t,{i_,j_,k_,l_},{x___}]] :=
        -Tensor[t,{j,i,k,l},{x}] /; !IndexOrderedQ[{i,j}];

  (** RULE 2
      The next rule specifies that

         i
      T        = 0  ,
       i   . .

      if "i" is an index associated with a Real bundle.  Note that
         _
         a            a
      T_        = -T         ,
       a   . .      a   . .

      if "a" is an index associated with a complex bundle by Rule 1 and
      the general rules for tensors.  **)

  t/: HoldPattern[Tensor[t,{L[b_],U[b_],_,_},{___}]] := 0 /; Type[Bundle[b]]===Real;

  (** RULE 3
      This rule says that

        T        = 0
         i i k l

      (This is relevant only for 1-dimensional bundles.)  **)

  t/: HoldPattern[Tensor[t,{i_,i_,_,_},{___}]] := 0;

  (** RULE 4
      This rule uses the symmetry

      T        = -T
       i j k l     i j l k

      to order the last two indices.  **)

  t/: HoldPattern[Tensor[t,{i_,j_,k_,l_},{x___}]] :=
        -Tensor[t,{i,j,l,k},{x}] /; !IndexOrderedQ[{k,l}];

  (** RULE 5
      This rule is the analogue of RULE 2 for the last two indices.  **)

  t/: HoldPattern[Tensor[t,{_,_,L[b_],U[b_]},{___}]] := 0 /; Type[Bundle[b]]===Real;

  (** RULE 6
      This is the analogue of RULE 3 for the last two indices.  **)

  t/: HoldPattern[Tensor[t,{_,_,i_,i_},{___}]] := 0;

  (** RULE 7
      This rule uses the symmetry

      T        = T
       i j k l    k l i j

      to order the first and third indices.  **)

  t/: HoldPattern[Tensor[t,{i_,j_,k_,l_},{x___}]] :=
        Tensor[t,{k,l,i,j},{x}] /; !IndexOrderedQ[{i,k}];

  (** RULE 8
      This rule takes care of the following problem.  Observe that

           i          i
      T        - T        = 0,
       i l   k    i k   l

      provided the bundle, for which i is an index, is real.  However,
      RULES 1-7 will not transform

           i                i
      T         into  T         .
       i l   k          i k   l

      This is why we have to define RULE 8.
      The result of RULE 8, together with the preceding rules, is that
      if an index in the first two is paired with an index in the last two,
      then the remaining indices will be ordered.  **)

  t/: HoldPattern[Tensor[t,{L[b_],j_,U[b_],l_},{x___}]]:=
        Tensor[t,{L[b],l,U[b],j},{x}] /;
        !IndexOrderedQ[{j,l}] && Type[Bundle[b]]===Real

  )

defSymmetryRelations[Tensor[t_,{},{}], TensorSymmetry[_,_,symlist_] ] :=
  Module[{i,len,perm,sgn,order,count,p},
    len = Length[symlist];
    For[i = 1, i < len, i++,
        perm = symlist[[i]];
        sgn  = symlist[[++i]];
        For[order = MyPermutationOrder[perm]; count = 1; p = perm,
            count < order,
            count++; p = p[[perm]],
            defTenSym[ Tensor[t,{},{}], p, sgn^count ];
       ];
    ];
  ];

defTenSym[Tensor[t_,{},{}],p_,sgn_] := (
   t/: HoldPattern[Tensor[t,{i__},{j___}]] :=
       sgn * Tensor[t,{i}[[p]],{j}] /;
       Length[{i}] == Length[p] && !IndexOrderedQ[ {i},{i}[[p]] ];
   If[Type[t]==={Complex},
      t/: HoldPattern[Tensor[t[Bar],{i__},{j___}]] :=
          sgn * Tensor[t[Bar],{i}[[p]],{j}] /;
          Length[{i}] == Length[p] && !IndexOrderedQ[ {i},{i}[[p]] ]]
   );

(****************** DefineTensorSymmetries *****************************)

SetAttributes[TensorSymmetry,HoldFirst];
SetAttributes[DefineTensorSymmetries,HoldFirst];

string[x_] := ToString[x];
string[x_,y__] := StringJoin[ ToString[x], ",", string[y] ];

Options[DefineTensorSymmetries] := {Quiet -> $Quiet};

DefineTensorSymmetries[name_, perms_, opts___] := Module[{quiet},
  If[ !optionsOK[DefineTensorSymmetries,opts], Return[] ];
  quiet = Quiet /. {opts} /. Options[DefineTensorSymmetries];
  (** Check that arguments are of the correct form **)
  Which[
   !symOK[perms],
     Message[DefineTensorSymmetries::invalid, HoldForm[perms],
                                              "tensor symmetry"],
   !newsymQ[name],
     Message[DefineTensorSymmetries::invalid, HoldForm[name],
                                              "tensor symmetry name"],
   True,
     If[!quiet, Print["Degree ", Length[perms[[1]]], " tensor symmetry ",
                       name," defined."]];
     name = TensorSymmetry[name,Length[perms[[1]]],perms];
   ];
 ];

symOK[ perms_List ] := permOK[ Length[ perms[[1]] ], perms ];

symOK[_] := False;

permOK[d_, {perm_List, sign_?ConstantQ, others___}] :=
  Length[perm]===d &&
  NontrivialPermutationQ[perm] &&
  permOK[ d, {others} ];

permOK[_, {} ] := True;

permOK[_, _] := False;

(****************** UndefineTensorSymmetries ****************************)

SetAttributes[UndefineTensorSymmetries,{HoldFirst,Listable}];

UndefineTensorSymmetries[name_] :=
  If[SameQ[ Head[name], TensorSymmetry ],
    Print["Tensor symmetry ",name," undefined."];
    ClearAll[name]
  ];

PermutationQ[{a__Integer}] := SameQ[ Range[Length[{a}]] , Sort[{a}] ];

NontrivialPermutationQ[{a__Integer}] :=
  Module[{ints},
    ints = Range[Length[{a}]];
    Return[ !SameQ[{a},ints] && SameQ[Sort[{a}],ints] ]
  ];

MyPermutationOrder[p_?PermutationQ] :=
  Module[{ints, perm, order},
    For[perm = p; ints = Range[Length[perm]]; order = 1,
        !SameQ[ints,perm],
        order++,
        perm = perm[[p]]
    ];
    Return[order]
  ];

(************************* Standard tensor definitions ********************)

DefineTensor[Curv,{2,2},Symmetries->{NoSymmetries,Skew},
  Quiet->True,Bundle->Automatic,Variance->{Co,Con,Co,Co}];

DefineTensor[Tor,{1,2},Symmetries->{NoSymmetries,Skew},
  Quiet->True,Bundle->Automatic,Variance->{Con,Co,Co}];

DefineTensor[Kronecker,2,Symmetries->Symmetric,Quiet->True,Bundle->Any,
                         Variance->{Covariant,Contravariant},
                         TeXFormat->"\\delta"];

DefineTensor[Basis,{1,1},Quiet->True,Bundle->Same,Variance->Same];

DefineTensor[Conn,{2,1},Quiet->True,Bundle->Automatic,Variance->{Co,Con,Co}];

Unprotect[Curv,Tor,Kronecker,Basis,Conn];

(* If the connection is flat, the curvature vanishes. *)

Curv/: HoldPattern[Tensor[Curv,{i_,j_,k___},{l___}]] := 0 /;
         FlatConnection[Bundle[i]];

(* The curvature matrix is skew-symmetric. *)

Curv/: HoldPattern[Tensor[Curv,{i_,j_,k___},{l___}]] :=
         - Tensor[Curv,{j,i,k},{l}] /;
         !IndexOrderedQ[{i,j}];

Curv/: HoldPattern[Tensor[Curv,{L[i_],U[i_],k___},{l___}]] := 0 /;
         Type[Bundle[i]] === Real;

(* The default connection is assumed compatible with each bundle.
   Thus Curv[L[i],U[k]] = Conn[L[i],U[k]] = 0 if i & k come from
   different bundles. *)

Curv/: HoldPattern[Tensor[Curv,{i_,j_,___},{___}]] := 0 /;
         Bundle[Lower[i]] =!= Conjugate[Bundle[Lower[j]]];

Conn/: HoldPattern[Tensor[Conn,{i_,j_,___},{___}]] := 0 /;
         Bundle[Lower[i]] =!= Conjugate[Bundle[Lower[j]]];

(* If the default frame is parallel, the connection coefficients vanish. *)

Conn/: HoldPattern[Tensor[Conn,{i_,j_,k___},{}]] := 0 /;
         ParallelFrame[Bundle[i]];

(* If the default frame is orthonormal, the connection forms are skew *)

Conn/: HoldPattern[Tensor[Conn,{i_,j_,k___},{}]] :=
         - Tensor[Conn,{j,i,k},{}] /;
         OrthonormalFrame[Bundle[i]] && !IndexOrderedQ[{i,j}];

(* If the default frame commutes, then the connection coefficients are
   symmetric modulo torsion. *)

Conn/: HoldPattern[Tensor[Conn,{i_,j_,k_},{}]] :=
         Tensor[Conn,{k,j,i},{}] + Tensor[Tor,{j,k,i},{}] /;
         Bundle[i]===Bundle[k] &&
         CommutingFrame[Bundle[i]] && !IndexOrderedQ[{i,k}];

(* "Covariant derivatives" of the connection coefficients never show up--
   convert them to ordinary directional derivatives. *)

Conn/: HoldPattern[Tensor[Conn,{i_,j_,k_},{l__}]] :=
  expandedCovD[ Tensor[Conn,{i,j,k},{}], {l} ];

(* If the connection is torsion-free, the torsion vanishes. *)

Tor/: HoldPattern[Tensor[Tor,{i_,j___},{k___}]] := 0 /;
         TorsionFree[Bundle[i]];

(* Rules for the Kronecker delta function (identity tensor) *)

(* Treat Kronecker as a metric when absorbing metrics in TensorSimplify *)

Kronecker/: HoldPattern[TensorMetricQ[Kronecker]] = True;

(* Kronecker is always parallel. *)

Kronecker/: HoldPattern[Tensor[Kronecker,{i_,j_},{k__}]] := 0;

(* Trace[Kronecker] = dimension *)

Kronecker/: HoldPattern[Tensor[Kronecker,{L[i_],U[i_]},{}]] :=
            Dimension[Bundle[i]];

(* Kronecker[L[i],L[j]] = metric[L[i],L[j]] if bundles match, 0 otherwise *)

Kronecker/: HoldPattern[Tensor[Kronecker,{LorU_[i_],LorU_[j_]},{}]] :=
            Metric[Bundle[i]][LorU[i],LorU[j]];

Kronecker/: HoldPattern[Tensor[Kronecker,{L[i_],U[j_]},{}]] := 0 /;
            !(Bundle[i] === Bundle[j]);

Kronecker/: HoldPattern[Tensor[Kronecker,{U[i_],L[j_]},{}]] := 0 /;
            !(Bundle[i] === Bundle[j]);

(* On a one-dimensional bundle, Kronecker[L[i],U[j]] = 1, even if i and
   j have different names *)

Kronecker/: HoldPattern[Tensor[Kronecker,{L[i_],U[j_]},{}]] := 1 /;
            Dimension[Bundle[i]] === 1;            

Kronecker/: HoldPattern[Tensor[Kronecker,{U[i_],L[j_]},{}]] := 1 /;
            Dimension[Bundle[i]] === 1;            

(* Rules for absorbing Kronecker in products.
   This is done more completely by AbsorbMetrics. *)

Tensor/: HoldPattern[Tensor[Kronecker,{a___,L[b_],c___},{}] *
         Tensor[t_,x___,{i___,U[b_],j___},y___]] :=
           Tensor[t,x,{i,a,c,j},y];

Tensor/: HoldPattern[Tensor[Kronecker,{a___,U[b_],c___},{}] *
         Tensor[t_,x___,{i___,L[b_],j___},y___]] :=
           Tensor[t,x,{i,a,c,j},y];

(* Basis[i] [j] gives Kronecker[i,j] (which will be changed to metric[i,j]
   if both indices have the same altitude. *)

Basis/: HoldPattern[Tensor[Basis,{i_,j_},{k___}]] :=
         Tensor[Kronecker,{i,j},{k}];

Protect[Curv,Tor,Kronecker,Basis,Conn];

(*************************** TensorExpr **********************)

(************************ Rank, Bundles, Variance *********************)

(* These three functions must know about EVERY function known to Ricci
   which might be held unevaluated as part of a tensor expression.  *)

(************************ Rank ****************************************)

(* There are two Rank functions:
   The public function "Rank" evaluates its argument, and then calls
   the Private function rank, which is a HoldAll function.  The reason
   is that rank is recursive, and since it's called so frequently,
   repeatedly re-evaluating its argument could cost too much time. *)

SetAttributes[ rank, HoldFirst ];

Rank[x_] := rank[x];

rank[HoldPattern[Tensor[Inverse[x_],{i_,j_},{}]]] := 0;
rank[HoldPattern[Tensor[t_[Bar],i_,j_]]] := rank[Tensor[t,i,j]];
rank[Tensor[t_,{},{}]] := TotalRank[t];
rank[Tensor[t_,{i___},{j__}]] = 0;
rank[Tensor[t_,{i__},{}]] := TotalRank[t] - Length[{i}];

rank[HoldPattern[Summation[x_]]] := rank[x];

rank[HoldPattern[Plus[x_,y__]]] := rank[x];

rank[HoldPattern[(x_)^(p_)]] := 0 /; 
  rank[x]===0 && rank[p]===0;
rank[HoldPattern[(x_)^(p_)]] := p rank[x] /; 
  IntegerQ[Unevaluated[p]] && p > 0;

rank[HoldPattern[Times[t__]]] := sumrank[t];
rank[HoldPattern[TensorProduct[t__]]] := sumrank[t];
rank[HoldPattern[Wedge[t__]]] := sumrank[t];

sumrank[t_,u__] := rank[t] + sumrank[u];
sumrank[t_] := rank[t];

rank[HoldPattern[Del[t_,opts___Rule]]] := rank[t] + 1;
rank[HoldPattern[Del[_,t_,opts___Rule]]] := rank[t];
rank[HoldPattern[Grad[t_,opts___Rule]]] := rank[t] + 1;
rank[HoldPattern[Div[t_,opts___Rule]]] := rank[t] -1;
rank[HoldPattern[Extd[t_]]] := rank[t] + 1;
rank[HoldPattern[ExtdStar[t_,opts___Rule]]] := rank[t] - 1;
rank[HoldPattern[Lie[_,t_]]] := rank[t];

rank[HoldPattern[Alt[t_]]] := rank[t];
rank[HoldPattern[Sym[t_]]] := rank[t];
rank[HoldPattern[Int[t_,u_]]] := rank[u]-rank[t];
rank[HoldPattern[Inner[t_,u_,opts___Rule]]] := 0;
rank[HoldPattern[HodgeInner[t_,u_,opts___Rule]]] := 0;
rank[HoldPattern[Dot[t___]]] := 
  Plus @@ (rank /@ {t}) - 2 (Length[Unevaluated[{t}]]-1);
rank[HoldPattern[Inverse[x_]]] := 2 /; rank[x]===2;
rank[HoldPattern[Det[x_]]] := 0 /; rank[x]===2;
rank[HoldPattern[Tr[x_]]] := 0 /; rank[x]===2;
rank[HoldPattern[Transpose[x_]]] := rank[x];
rank[HoldPattern[Conjugate[x_]]] := rank[x];

rank[_RiemannTensor] := 4;
rank[_RicciTensor]   := 2;
rank[_ScalarCurv]    := 0;
rank[_LeviCivitaConnection] := 3;
rank[_Curvature]     := 4;

rank[f_[g_]] := 0 /; mathFunctionQ[f] && rank[g] === 0;

rank[x_] := 0 /; ConstantQ[Unevaluated[x]];

(* The following rule is for TimesFormat, which needs to compute
   the rank of an expression wrapped in HoldForm *)

rank[ HoldPattern[ HoldForm[x_] ] ] := rank[x];

(*********** Functions for determining the type of a tensor expression *****)

ScalarQ[x_] := 
  rank[x]===0 && GetFreeIndices[Unevaluated[x]]==={};

VectorFieldQ[x_] := Variance[Unevaluated[x]]==={Contravariant};
                    (* This will be true only if rank[x] = 1 *)

FormQ[x_] := 
  SkewQ[Unevaluated[x]] && FreeQ[Variance[Unevaluated[x]],Contravariant];

(************************* SymmetricQ *******************************)

SymmetricQ[Tensor[t_[Bar],i_,j_]]    := 
  SymmetricQ[ Conjugate[ Tensor[t[Bar],i,j] ] ];
SymmetricQ[Tensor[t_,i_,j_][args__]] := SymmetricQ[Tensor[t,i,j]];
SymmetricQ[Tensor[t_,{i___},{j___}]] := True /; 
  rank[Tensor[t,{i},{j}]] <= 1;
SymmetricQ[Tensor[t_,{},{}]]         := 
  MemberQ[ {Symmetric,Hermitian}, Symmetries /. TensorData[t]];

(** If the tensor is partially indexed, call it symmetric if it
    has two symmetries and the second one is Symmetric. **)

SymmetricQ[Tensor[t_,{i__},{}]] := 
  MatchQ[Symmetries /. TensorData[t], {_,Symmetric}|{_,Hermitian}];

SymmetricQ[HoldPattern[Times[x__]]]  := And @@ SymmetricQ /@ {x};
SymmetricQ[HoldPattern[Plus[x__]]]   := And @@ SymmetricQ /@ {x};
SymmetricQ[x_ ^ p_]  := (SymmetricQ[x] && IntegerQ[p] && p > 0) ||
                        (rank[x]===0 && rank[p]===0);
SymmetricQ[HoldPattern[Summation[x_]]] := SymmetricQ[x];
SymmetricQ[_Sym]     := True;
SymmetricQ[HoldPattern[Div[x_,opts___Rule]]] := SymmetricQ[x];
SymmetricQ[HoldPattern[Del[_,x_,opts___Rule]]] := SymmetricQ[x];
SymmetricQ[HoldPattern[Lie[_,x_]]] := SymmetricQ[x];
SymmetricQ[HoldPattern[Derivative[k__][t_Tensor][args__]]] :=
  SymmetricQ[t];

SymmetricQ[HoldPattern[Dot[x_,y_]]] := SymmetricQ[x] /; rank[y]===1 && rank[x]>2;
SymmetricQ[HoldPattern[Dot[x_,y_]]] := SymmetricQ[y] /; rank[x]===1 && rank[y]>2;
SymmetricQ[HoldPattern[Inverse[x_]]] := SymmetricQ[x];
SymmetricQ[HoldPattern[Conjugate[x]]] := SymmetricQ[x];

SymmetricQ[HoldPattern[_RicciTensor]] := True;

SymmetricQ[x_]       := (rank[x]===0 || rank[x]===1);

(******************** SkewQ[x] is True iff x is skew-symmetric *******)

SkewQ[Tensor[t_[Bar],i_,j_]]    := 
  SkewQ[ Conjugate[ Tensor[t[Bar],i,j] ] ];
SkewQ[Tensor[t_,i_,j_][args__]] := SkewQ[Tensor[t,i,j]];
SkewQ[Tensor[t_,{i___},{j___}]] := True /; 
  rank[Tensor[t,{i},{j}]] <= 1;
SkewQ[Tensor[t_,{},{}]] := 
  MemberQ[{Alternating,SkewHermitian}, Symmetries /. TensorData[t]];

(** If the tensor is partially indexed, call it skew if it
    has two symmetries and the second one is Alternating. **)

SkewQ[Tensor[t_,{i__},{}]] := 
  MatchQ[Symmetries /. TensorData[t], {_,Alternating}|{_,SkewHermitian}];

SkewQ[HoldPattern[Plus[x__]]]  := And @@ SkewQ /@ {x};
SkewQ[f_ * g_] := SkewQ[g] /; rank[f]===0;
SkewQ[x_ ^ p_] := (rank[x]===0 && rank[p]===0);
SkewQ[_Wedge] := True;
SkewQ[_Alt]   := True;
SkewQ[_Extd]  := True;
SkewQ[_ExtdStar] := True;
SkewQ[HoldPattern[Div[x_,opts___Rule]]] := SkewQ[x];
SkewQ[HoldPattern[Del[_,x_,opts___Rule]]] := SkewQ[x];
SkewQ[HoldPattern[Lie[_,x_]]] := SkewQ[x];
SkewQ[HoldPattern[Derivative[k__][t_Tensor][args__]]] :=
  SkewQ[t];

SkewQ[HoldPattern[_Int]] := True;
SkewQ[HoldPattern[Dot[x_,y_]]] := SkewQ[x] /; rank[y]===1 && rank[x]>2;
SkewQ[HoldPattern[Dot[x_,y_]]] := SkewQ[y] /; rank[x]===1 && rank[y]>2;
SkewQ[HoldPattern[Inverse[x_]]] := SkewQ[x];
SkewQ[HoldPattern[Conjugate[x]]] := SkewQ[x];

SkewQ[x_]     := (rank[x]===0 || rank[x]===1);

(*********************** Variance ***********************************)

(* Variance[exp] returns a list whose length is equal to Rank[exp], and
   each of whose entries is Covariant or Contravariant. *)

Unprotect[Variance];
HoldPattern[Variance[Tensor[Inverse[x_],{i_,j_},{}]]] := {};
Variance[Tensor[t_,{},{}]] := Variance /. TensorData[t] //. {Same->Covariant};
Variance[Tensor[t_,{i__},{}]] := 
  Drop[Variance /. TensorData[t], Length[{i}]] //.
  {Same -> Switch[IndexAltitude[ Last[{i}] ],
                  Upper, Covariant,
                  Lower, Contravariant]};
Variance[Tensor[t_[Bar],i_,j_]] := Variance[ Conjugate[ Tensor[t[Bar],i,j] ] ];
Variance[Tensor[t_,i_,j_][args__]] := Variance[Tensor[t,i,j]];

Variance[HoldPattern[Plus[x__]]] := Variance[First[{x}]];
Variance[HoldPattern[Summation[x_]]] := Variance[x];
Variance[HoldPattern[x_^q_]] := 
  Flatten[Table[Variance[x],{q}]] /; IntegerQ[q] && q > 0;
Variance[HoldPattern[TensorProduct[t__]]] := 
  Flatten[ Variance[#]& /@ {t}];
Variance[HoldPattern[Wedge[t__]]] := 
  Flatten[ Variance[#]& /@ {t}];
Variance[HoldPattern[Times[t__]]] := 
  Flatten[ Variance[#]& /@ {t}];
Variance[HoldPattern[Del[t_,opts___Rule]]] :=  
  Append[ Variance[t], Covariant ];
Variance[HoldPattern[Grad[t_,opts___Rule]]] :=  
  Append[ Variance[t], Contravariant ];
Variance[HoldPattern[Del[_,t_,opts___Rule]]] := Variance[t];
Variance[HoldPattern[Lie[_,t_]]] := Variance[t];
Variance[HoldPattern[Div[t_,opts___Rule]]] := Drop[Variance[t],-1];
Variance[HoldPattern[Extd[t_]]] := 
  Append[ Variance[t], Covariant];
Variance[HoldPattern[ExtdStar[t_,opts___Rule]]] := Drop[Variance[t],-1];
Variance[HoldPattern[Derivative[k__][t_Tensor][args__]]] :=
  Variance[t];

Variance[HoldPattern[Alt[x_]]] := Variance[x];
Variance[HoldPattern[Sym[x_]]] := Variance[x];
Variance[HoldPattern[Int[x_,y_]]] := Drop[Variance[y],rank[x]];
Variance[HoldPattern[Dot[x_,y__]]] := 
  Join[ Drop[Variance[x],-1], Drop[Variance[Dot[y]],1] ];
Variance[_U] := Contravariant;
Variance[_L] := Covariant;
Variance[x_] := {} /; rank[x]===0;
Variance[HoldPattern[Inverse[x_]]] :=
  Pair /@ Reverse[ Variance[x] ];
Variance[HoldPattern[Transpose[x_]]] :=
  Reverse[ Variance[x] ];
Variance[HoldPattern[Conjugate[x_]]] := Variance[x];

Variance[_RiemannTensor] := Table[Covariant,{4}];
Variance[_RicciTensor] := Table[Covariant,{2}];
Variance[_LeviCivitaConnection] := {Covariant,Contravariant,Covariant};
Variance[_Curvature] := {Co,Con,Co,Co};

Pair[Contravariant] := Covariant;
Pair[Covariant]     := Contravariant;

Protect[Variance];

(************************** Bundles ******************************)

(* Bundles[exp] returns a list whose length is equal to Rank[x], and
   each of whose entries is a sublist, denoting the allowable bundles
   for that index slot.  Special values: {Any} means the slot will
   accept an index from any bundle; {Null} means no bundle was 
   specified. *)

Bundles[Tensor[Inverse[x_],{i_,j_},{}]] := {};
Bundles[Tensor[t_[Bar],i_,j_]]  := 
  Conjugate[Bundles[Conjugate[Tensor[t[Bar],i,j]]]];

Bundles[Tensor[t_,{},{}]] :=   
  Bundles /. TensorData[t] //. {Automatic -> Any, Same -> Any};

Bundles[Tensor[t_,{i__},{}]] := 
  (Drop[Bundles /. TensorData[t], Length[{i}]] //. 
    {{Automatic} -> TangentBundle[ Bundle[ First[{i}] ] ],
     {Same} -> Union[{#},Conjugate[{#}]]& @ Bundle[ Last[{i}] ] });  

Bundles[Tensor[t_,{i___},{j__}]] := {};

Bundles[HoldPattern[ x_ + Tensor[Conn,{},{}] ]] := Bundles[x];
Bundles[HoldPattern[x_Plus]] := 
  (Union @@ # &) /@ Transpose[ Bundles /@ List @@ x ] //.
  { {a___,Any,b___} :> {a,b} /; Length[{a,b}] > 0 }; 
Bundles[HoldPattern[Summation[x_]]] := Bundles[x];
Bundles[HoldPattern[x_^q_]] := 
  Flatten[Table[Bundles[x],{q}], 1] /; IntegerQ[q] && q > 0;
Bundles[HoldPattern[TensorProduct[t__]]] := 
  Flatten[ Bundles /@ {t}, 1 ];
Bundles[HoldPattern[Wedge[t__]]] := 
  Table[ Union[ Flatten[ Bundles /@ {t} ]],
         {rank[Wedge[t]]}
       ];
Bundles[HoldPattern[Times[t__]]] := 
  Table[ Union[ Flatten[ Bundles /@ {t} ]],
         {rank[Wedge[t]]}
       ];

Bundles[HoldPattern[Del[t_,con__Rule]]] :=  
  Append[ Union[connectionbundle[con],#]& /@ Bundles[t], 
          UnderlyingTangentBundle[t] ];
Bundles[HoldPattern[Grad[t_,con__Rule]]] :=  
  Append[ Union[connectionbundle[con],#]& /@ Bundles[t], 
          UnderlyingTangentBundle[t] ];
Bundles[HoldPattern[Del[_,t_,con__Rule]]] := 
  Union[connectionbundle[con],#]& /@ Bundles[t];
Bundles[HoldPattern[Div[t_,con__Rule]]] := 
  Drop[Union[connectionbundle[con],#]& /@ Bundles[t],-1];
Bundles[HoldPattern[ExtdStar[t_,con__Rule]]] := 
  Drop[Union[connectionbundle[con],#]& /@ Bundles[t],-1];

Bundles[HoldPattern[Del[t_]]] :=  
  Append[ Bundles[t], UnderlyingTangentBundle[t] ];
Bundles[HoldPattern[Grad[t_]]] :=  
  Append[ Bundles[t], UnderlyingTangentBundle[t] ];
Bundles[HoldPattern[Del[_,t_]]] := Bundles[t];
Bundles[HoldPattern[Lie[_,t_]]] := 
  Table[ UnderlyingTangentBundle[t], {rank[t]} ];
Bundles[HoldPattern[Div[t_]]] := Drop[Bundles[t],-1];
Bundles[HoldPattern[Extd[t_]]] := 
  Table[ UnderlyingTangentBundle[t], {rank[t]+1} ];
Bundles[HoldPattern[ExtdStar[t_]]] := Drop[Bundles[t],-1];

connectionbundle[Connection -> HoldPattern[Tensor[Conn,{},{}] + diff_]] :=
  Bundles[diff] [[2]];
connectionbundle[Connection -> con_] := 
  Bundles[con] [[2]];
connectionbundle[Metric -> met_] := Bundles[met][[1]];
connectionbundle[__, Connection -> con_, ___] :=
  connectionbundle[Connection -> con];
connectionbundle[Connection -> con_, __] :=
  connectionbundle[Connection -> con];

Bundles[HoldPattern[Alt[x_]]] := Table[ Union @@ Bundles[x], {rank[x]} ];
Bundles[HoldPattern[Sym[x_]]] := Table[ Union @@ Bundles[x], {rank[x]} ];
Bundles[HoldPattern[Int[x_,y_]]] := Drop[Bundles[y],rank[x]];
Bundles[HoldPattern[Dot[x_,y__]]] := 
  Join[ Drop[Bundles[x],-1], Drop[Bundles[Dot[y]],1] ];
Bundles[HoldPattern[Inverse[x_]]] := Reverse[Bundles[x]];
Bundles[HoldPattern[Transpose[x_]]] := Reverse[Bundles[x]];
Bundles[HoldPattern[Conjugate[x_]]] := (Conjugate /@ # &) /@ Bundles[x];

Bundles[HoldPattern[RiemannTensor[h_]]] := Table[Bundles[h][[1]],{4}];
Bundles[HoldPattern[RicciTensor[h_]]]   := Bundles[h];
Bundles[HoldPattern[LeviCivitaConnection[h_]]] := Table[Bundles[h][[1]],{3}];
Bundles[HoldPattern[Curvature[ c_ + Tensor[Conn,{},{}] ]]] := 
  Join[ Table[ Bundles[c] [[1]], {2} ], 
        Table[ UnderlyingTangentBundle[c], {2} ] ];
Bundles[HoldPattern[Curvature[c_]]] := 
  Join[ Table[ Bundles[c] [[1]], {2} ], 
        Table[ UnderlyingTangentBundle[c], {2} ] ];

Bundles[x_] := {} /; rank[x]===0;

(****************** UnderlyingTangentBundle *************************)

(* This function accepts any tensor expression, and returns a list
   of bundles, representing the direct-sum decomposition of the
   underlying tangent bundle of the expression. *)

UnderlyingTangentBundle[Tensor[t_,{},{}]] := 
  TangentBundle /. TensorData[t] /. {Null} -> Flatten[{$DefaultTangentBundle}];

UnderlyingTangentBundle[x_] := Which[
  ConstantQ[x],        Flatten[{$DefaultTangentBundle}],
  !(FreeQ[x,L[_]]) || !(FreeQ[x,U[_]]),
                       TangentBundle @ Bundle @ First @ GetIndices @ x,
  rank[x] > 0,         TangentBundle[ Bundles[x] [[1,1]] ],
  True,                UnderlyingTangentBundle @ First @ Occurrences[x,_Tensor]
  ];

(*************************** Sym **********************************)

Sym[x_Plus] := Sym /@ x;
Sym[(f_ /; rank[f]===0) * t_] := f*Sym[t];
Sym[ HoldPattern[TensorProduct[x__]] ] := Times @@ Sym /@ {x};
Sym[ HoldPattern[Transpose[x_]] ] := Sym[x];

Sym[x_?SymmetricQ] := x;
Sym[x_?SkewQ]      := 0;

Sym/: TCompute[ HoldPattern[Sym[f_]], {ind___} ] :=
  Module[{perms},
    perms = Permutations[{ind}];
    Return[Plus @@ (InsertIndices[ f, # ]& /@ perms)  / Length[perms] ]
    ];

Sym/: Conjugate[HoldPattern[Sym[x_]]] := Sym[Conjugate[x]];

(************************** Alt **************************************)

Alt[x_Plus] := Alt /@ x;
Alt[(f_ /; rank[f]===0) * t_] := f*Alt[t];
Alt[ HoldPattern[TensorProduct[x__]] ] := 
  (1/WedgeFactor[rank/@{x}]) Wedge @@ Alt /@ {x};
Alt[ HoldPattern[Transpose[x_]] ] := (-1)^(rank[x](rank[x]-1)/2) * Alt[x];

Alt[x_?SkewQ]      := x;
Alt[x_?SymmetricQ] := 0;

(** Here's where we insert indices **)

Alt/: TCompute[ HoldPattern[Alt[f_]], {ind__} ] :=  
  Module[{perms},
    perms = Permutations[{ind}];
    Return[Plus @@ (Signature[{ind}] * (Signature /@ perms)
            * (InsertIndices[ f, # ]& /@ perms) ) / Length[perms] ]
    ];

Alt/:       Conjugate[HoldPattern[Alt[x_]]] := Alt[Conjugate[x]];

(****************** Power *********************************)

Unprotect[Power];

(* The Summation function is inserted to prevent (a b)^p from being
   expanded into a^p b^p when p is an integer and a,b contain indices. *)

Power/: HoldPattern[Times[a__]^p_] := 
  Times @@ (# ^ p&) /@ Select[ {a},  FreeQ[#,L|U]& ] *
  Summation[  Times @@ Select[ {a}, !FreeQ[#,L|U]& ] ]^p;

Power/: TCompute[ HoldPattern[Power[ x_, p_ ]], {} ] :=
  InsertIndices[ x, {} ] ^ InsertIndices[ p, {} ];

Power/: TCompute[ HoldPattern[Power[x_,p_]], {i___} ] := Module[ {n},
  Which[
    !SymmetricQ[x],
       Print["ERROR: non-symmetric tensor used in symmetric product"];
       ERROR[InsertIndices[x^p,{i}]],
    rank[x] > 0 && GetDummyIndices[x] =!= {},
       InsertIndices[ Times @@ Table[ NewDummy[x], {n,p} ], {i} ],
    True,
       (* The argument must be Unevaluated below, because otherwise
          Sym[a(X)a] would be converted back to a*a *)
       TCompute[ Unevaluated[ Sym[TensorProduct @@ Table[x,{n,p}]]], {i} ]
    ]
  ];

Protect[Power];

(*************************** Summation ******************************)

(* Mathematica automatically expands (a b)^p into a^p b^p when p is an 
   integer.  When a and/or b are tensor expressions with indices, this
   is not what we want.  In such cases, (a b)^p is converted to 
   Summation[a b]^p to inhibit the expansion. (See Power above.)

   The following rules are meant to get rid of the Summation function 
   if it ever shows up in an algebraic expression anywhere except 
   raised to a power. *)
 
Summation/: HoldPattern[Summation[x_]^1] := x;
Summation/: HoldPattern[Summation[x_] + y_] := x+y;
Summation/: HoldPattern[Summation[x_] * y_] := x*y;
HoldPattern[Summation[x_ /; (Head[x] =!= Times || FreeQ[x, (L|U)]) ]] := x;
HoldPattern[Summation[(a_ /; FreeQ[a, (L|U)]) * b_]] := a*Summation[b];
Summation/: Conjugate[ HoldPattern[Summation[x_]] ] := Summation[Conjugate[x]];

Summation/: TCompute[ HoldPattern[Summation[x_]], {i___} ] :=
  TCompute[ x, {i} ];

(******************** Inverse ***********************************)

Unprotect[Inverse];

HoldPattern[Inverse[Inverse[x_]]] := x;
HoldPattern[Inverse[f_ x_]] := (1/f) Inverse[x] /; Rank[f]===0 && !PairedQ[f,x];
HoldPattern[Inverse[Dot[x___]]] := Dot @@ Reverse[ Inverse /@ {x} ];
HoldPattern[Inverse[Tensor[Kronecker,{},{}]]] := Kronecker;
HoldPattern[Inverse[a_ + b_]] := Inverse[a] + Inverse[b] /;
  Intersection[ Union @@ Bundles[a], Union @@ Bundles[b] ] === {};
Inverse/: HoldPattern[Dot[ a___, Inverse[x_], x_, b___ ]] := Dot[a,Kronecker,b];
Inverse/: HoldPattern[Dot[ a___, x_, Inverse[x_], b___ ]] := Dot[a,Kronecker,b];

Inverse/: Conjugate[HoldPattern[Inverse[x_]]] := Inverse[Conjugate[x]];

(* If x has changed since we last looked at it, re-evaluate it and
   re-insert indices.  This is necessary since Tensor has the "HoldFirst"
   attribute. *)

Inverse/: Tensor[Inverse[x_], {i_,j_}, {}] := 
  Inverse[x] [i,j] /; ValueQ[Inverse[x]];

(* The indexed version of the inverse of a metric is just the metric
   itself *)

Inverse/: HoldPattern[Tensor[Inverse[g_], {i_,j_}, {}]] := 
  InsertIndices[g, {i,j}] /; MetricQ[g];

(* The symmetries of the inverse are those of the tensor itself. *)

Inverse/: HoldPattern[Tensor[Inverse[x_],{i_,j_},{}]] :=
  Tensor[Inverse[x],{j,i},{}] /; !IndexOrderedQ[{i,j}] && SymmetricQ[x];

Inverse/: HoldPattern[Tensor[Inverse[x_],{i_,j_},{}]] :=
  - Tensor[Inverse[x],{j,i},{}] /; !IndexOrderedQ[{i,j}] && SkewQ[x];

Inverse/: HoldPattern[Tensor[Inverse[Tensor[h_Symbol,{},{}]],{i_,j_},{}]] := 
  0/; (MemberQ[{Hermitian,SkewHermitian}, Symmetries /. TensorData[h]] &&
       SameQ @@ Bundle /@ Lower /@ {i,j});  

(* Here's where we compute the derivative of an indexed inverse tensor *)

Inverse/: HoldPattern[Tensor[Inverse[x_],{i_,j_},{k_,l___}]] :=
  - CovD[ InsertIndices[ Dot[ NewDummy[Inverse[x]],
                              Del[Basis[k],x],
                              NewDummy[Inverse[x]] ],
                         {i,j}
                       ],
          {l}];                            

Inverse/: TCompute[ HoldPattern[Inverse[x_]], {i__} ] := (
  Print["ERROR: ",x," is not a 2-tensor"];
  Return[ERROR[Inverse[x][i]]]
  ) /; !Rank[x]===2;

Inverse/: TCompute[ HoldPattern[Inverse[x_]], {i__} ] := 
  Tensor[Inverse[x],{i},{}];

Protect[Inverse];

(*********************************** Det *******************************)

Unprotect[Det];

Det[ HoldPattern[Transpose[x_]] ] := Det[x];
Det[f_ x_] := f ^ (Plus @@ (Dimension /@ First[Bundles[x]])) Det[x] /;
  rank[f]===0 && !MemberQ[ First[Bundles[x]], Any|Null];
Det[ HoldPattern[Dot[a___]] ] := Times @@ (Det /@ {a}) /;
  MatchQ[ rank /@ {a}, {(2)..} ];
Det[ HoldPattern[Inverse[h_]] ] := 1/Det[h];
Det[ HoldPattern[Tensor[g_,{},{}]] ] := 1 /;
  MetricQ[g];
Det/: Conjugate[ HoldPattern[Det[x_]] ] := Det[Conjugate[x]];

TCompute[ HoldPattern[Det[x_]], {}] := Det[x];

Protect[Det];

(***************************** Transpose *******************************)

Unprotect[Transpose];

MakeLinear[Transpose];
Transpose[ HoldPattern[Transpose[x_]] ] := x;
Transpose[f_ * x_] := f Transpose[x] /; rank[f]===0;
Transpose[ HoldPattern[ TensorProduct[x__] ] ] := 
  Reverse[Transpose /@ TensorProduct[x]];
Transpose[ HoldPattern[ Dot[x__] ] ] :=
  Reverse[Transpose /@ Dot[x]];
Transpose[x_] := x /; Head[x]=!=List && (SymmetricQ[x] || rank[x] < 2);
Transpose[x_] := (-1)^(rank[x](rank[x]-1)/2) * x /; Head[x]=!=List && SkewQ[x];
Transpose/: Conjugate[ HoldPattern[Transpose[x]] ] := Transpose[Conjugate[x]];

Transpose/: TCompute[ HoldPattern[Transpose[x_]], {i___}] := 
  TCompute[x, Reverse[{i}]];

Protect[Transpose];

(***************************** Tr **************************************)

Unprotect[Tr];
MakeLinear[Tr];
Tr[ Tensor[g_,{},{}] ] := Plus @@ (Dimension /@ First[Bundles[g]]) /;
  MetricQ[g];
Tr[ f_ x_ ] := f Tr[x] /; rank[f]===0;
Tr[ HoldPattern[Transpose[x_]] ] := Tr[x];
Tr/: Conjugate[ HoldPattern[Tr[x_]] ] := Tr[Conjugate[x]];
Tr[ HoldPattern[Del[v_,x_,con___Rule]] ] := Del[v, Tr[x], con];
Tr[ x_ ] := 0 /; SkewQ[x] && Rank[x]===2;

TCompute[ HoldPattern[Tr[x_]], {} ] := 
  Plus @@ (InsertIndices[ x, {L[#],U[#]} ]& /@ NewBundleSymbol /@ 
                                               First[Bundles[x]]);
Protect[Tr];

(************************ Math Functions *************************)

SetAttributes[DefineMathFunction,{Listable}];

Options[DefineMathFunction] := {Type -> Real,
                                 Quiet -> $Quiet};

DefineMathFunction[f_, opts___] :=
  Module[{type,quiet,protected},
    If[!optionsOK[DefineMathFunction,opts],Return[]];
    {type,quiet} = {Type,Quiet} /. {opts} /. Options[DefineMathFunction];
    Which[
      !funtypeQ[type ],
          Message[DefineMathFunction::invalid,type,"math function type"],
      Head[f] =!= Symbol,
          Message[DefineMathFunction::invalid,f,"math function"],
      True,
          protected = Unprotect[f];
          f/: HoldPattern[f[x_][i:(_L|_U)...]] := InsertIndices[f[x],{i}];
          $MathFunctions = Union[ $MathFunctions, {f} ];
          Protect[protected];
          DeclareMathFunction[f,Type->type];
          If[ !quiet, Print["Math function ",f," defined."] ];
      ];
    Return[f];
  ];

(* funtypeQ tests the Type option for DefineMathFunction and
   DeclareMathFunction to make sure it's valid and not
   self-contradictory *)

funtypeQ[t_List] := 
  (t === {Complex}) ||
  (t === {Imaginary}) ||
  (t === {Automatic}) ||
  (ContainsQ[ {Real,Positive,Negative,NonPositive,NonNegative}, t ] &&
   !ContainsQ[ t, {Positive,Negative} ] &&
   !ContainsQ[ t, {Positive,NonPositive} ] &&
   !ContainsQ[ t, {Negative,NonNegative} ]);
funtypeQ[t_] := funtypeQ[{t}];

Options[DeclareMathFunction] := {Type -> Null};

DeclareMathFunction[ f_, opts___ ] := 
  Module[{type,protected,vals,attrs,assignments},
  protected = Unprotect[f];
  If[MemberQ[optionNames[{opts}], Type], 
    type = Flatten[{Type /. {opts}}];
    If[ !funtypeQ[type], 
        Message[Declare::invalid, type, "math function type"];
        Return[] ];
    vals = First /@ UpValues[f];
    attrs = Position[ vals /. HoldPattern->Hold,
                         Hold[Sign[f[_]]] |
                         Hold[NonPositive[f[_]]] |
                         Hold[Negative[f[_]]] |
                         Hold[Positive[f[_]]] |
                         Hold[NonNegative[f[_]]] |
                         Hold[Conjugate[f[_]]] ];
    assignments = PartAt[vals,#]& /@ attrs;
    (TagUnset[f,#]&) /@ assignments;
    If[ MemberQ[type, Real|Positive|Negative|
                      NonNegative|NonPositive], 
        f/: HoldPattern[Conjugate[f[x_]]]   := f[x] ];
    If[ MemberQ[type, Imaginary],               
        f/: HoldPattern[Conjugate[f[x_]]]   := -f[x] ];
    If[ MemberQ[type, Positive],                
        f/: HoldPattern[Sign[f[x_]]]        := 1;
        f/: HoldPattern[Positive[f[x_]]]    := True;
        f/: HoldPattern[NonNegative[f[x_]]] := True  ];
    If[ MemberQ[type, NonNegative],    
        f/: HoldPattern[NonNegative[f[x_]]] := True;
        f/: HoldPattern[Negative[f[x_]]]    := False ];
    If[ MemberQ[type, NonPositive],    
        f/: HoldPattern[Positive[f[x_]]]    := False;
        f/: HoldPattern[NonPositive[f[x_]]] := True  ];
    If[ MemberQ[type, Negative],                
        f/: HoldPattern[Sign[f[x_]]]        := -1;
        f/: HoldPattern[Positive[f[x_]]]    := False;
        f/: HoldPattern[NonNegative[f[x_]]] := False ];
    If[ MemberQ[type, Automatic],               
        f/: HoldPattern[Conjugate[f[x_]]]   := f[Conjugate[x]] ];
    ];
  Protect[protected];
  Return[f]
  ];

mathFunctionQ[ f_Symbol ] := MemberQ[ $MathFunctions, f];
mathFunctionQ[ Derivative[i_][ f_Symbol ] ] := mathFunctionQ[f];
mathFunctionQ[ _ ] := False;

TCompute[ f_?mathFunctionQ[g_], {} ] := f[ InsertIndices[g,{}] ];

TCompute[ Conjugate[f_?mathFunctionQ[g_]], {} ] := Conjugate[f[ InsertIndices[g,{}] ]];

(***************** Conjugate, Re, and Im *******************************)

Unprotect[{Re,Im}];
Re[a_Plus] := Re /@ a;
Im[a_Plus] := Im /@ a;
Re[a_] := a/2 + Conjugate[a]/2;
Im[a_] := -I*a/2 + I*Conjugate[a]/2;
Protect[{Re,Im}];

(** We need to teach Conjugate a few more tricks **)

Unprotect[Conjugate];

Conjugate[Conjugate[x_]] := x;
Conjugate[x_Plus]  := Conjugate /@ x;
Conjugate[x_Times] := Conjugate /@ x;
HoldPattern[ Conjugate[x_ ^ (y_?IntegerQ)] ] := Conjugate[x]^y;
HoldPattern[ Conjugate[(x_?NonNegative) ^ y_] ] := x^Conjugate[y];
HoldPattern[ Conjugate[Log[x_?NonNegative]] ] := Log[x];
HoldPattern[ Conjugate[f_[g_]][]] := Conjugate[f[g]] /;
  mathFunctionQ[f];
HoldPattern[ Conjugate[f_[g_]][x__?IndexQ]] := CovD[Conjugate[f[g]],{x}] /;
  mathFunctionQ[f];

Protect[Conjugate];

(************************* Positive *********************************)

Unprotect[Positive];

HoldPattern[ Positive[ Tensor[t_Symbol,{},{}] /;
                   MemberQ[Type[t], Positive] ] ] := True;
HoldPattern[ Positive[ Tensor[t_Symbol,{},{}] /;
                   MemberQ[Type[t], NonPositive|Negative] ] ] := False;

HoldPattern[ Positive[ a_Plus /; (And @@ (NonNegative /@ (List @@ a))) &&
                             MemberQ[ Positive /@ (List @@ a), True ] ]] := True;
HoldPattern[ Positive[ a_Plus /; (And @@ (NonPositive /@ (List @@ a)))]] := False;
HoldPattern[ Positive[ (a_?Positive) b_ ] ] := Positive[b];
HoldPattern[ Positive[ (a_?Negative) b_ ] ] := Negative[b];

HoldPattern[ Positive[ (a_?Positive) ^ (b_ /; Conjugate[b]===b) ] ] := True;
HoldPattern[ Positive[ (a_?Negative) ^ (b_?EvenQ) ]] := True;
HoldPattern[ Positive[ (a_?NonPositive) ^ (b_?OddQ) ] ] := False;

Protect[Positive];

(************************* Sign *********************************)

Unprotect[Sign];

HoldPattern[ Sign[ Tensor[t_Symbol,{},{}] /;
                   MemberQ[Type[t], Positive] ] ] := 1;
HoldPattern[ Sign[ Tensor[t_Symbol,{},{}] /;
                   MemberQ[Type[t], Negative] ] ] := -1;

Protect[Sign];
(******************* NonNegative ********************************)

Unprotect[NonNegative];

HoldPattern[ NonNegative[ Tensor[t_Symbol,{},{}] /;
                      MemberQ[Type[t], 
                              NonNegative|Positive] ] ] := True;
HoldPattern[ NonNegative[ Tensor[t_Symbol,{},{}] /;
                      MemberQ[Type[t], 
                              Negative] ] ] := False;

HoldPattern[ NonNegative[ a_Plus /; (And @@ (NonNegative /@ (List @@ a))) ]] := True;
HoldPattern[ NonNegative[ a_Plus /; (And @@ (NonPositive /@ (List @@ a))) &&
                             MemberQ[ Negative /@ (List @@ a), True ] ]] := False;
HoldPattern[ NonNegative[ (a_?Positive) b_ ] ] := NonNegative[b];
HoldPattern[ NonNegative[ (a_?Negative) b_ ] ] := NonPositive[b];
HoldPattern[ NonNegative[ a_Times /; (And @@ (NonNegative /@ (List @@ a))) ]] := True;
HoldPattern[ NonNegative[ (a_ ^ p_.) * (Conjugate[a_] ^ p_.) * b_. 
                      /; Conjugate[p]===p ]] := True;

HoldPattern[ NonNegative[ (a_?NonNegative) ^ (b_ /; Conjugate[b]===b) ] ] := True;
HoldPattern[ NonNegative[ (a_ /; Conjugate[a]===a) ^ (b_?EvenQ) ] ] := True;
HoldPattern[ NonNegative[ (a_?Negative) ^ (b_?OddQ) ] ] := False;

HoldPattern[ NonNegative[ Tensor[t_, i_, j_] * Tensor[t_, k_, l_] * b_. /;
                   t===Conjugate[t] && 
                   Pair /@ k === i && 
                   Pair /@ l === j &&
                   And @@ PositiveDefinite /@ Bundle /@ Join[i,j]
                 ]]:= NonNegative[b];
HoldPattern[ NonNegative[ Tensor[t_, i_, j_] * Tensor[t_[Bar], k_, l_] * b_. /;
                   Pair /@ k === i && 
                   Pair /@ l === j &&
                   And @@ PositiveDefinite /@ Bundle /@ Join[i,j]
                 ]]:= NonNegative[b];

HoldPattern[ NonNegative[ Inner[ a_, b_ ] /; 
                   Conjugate[a]===b &&
                   And @@ PositiveDefinite /@ Flatten[Join[Bundles[a],Bundles[b]]]
                 ] ] := True;

HoldPattern[ NonNegative[ HodgeInner[ a_, b_ ] /; 
                   Conjugate[a]===b &&
                   And @@ PositiveDefinite /@ Flatten[Join[Bundles[a],Bundles[b]]]
                 ] ] := True;

Protect[NonNegative];

(******************* NonPositive ********************************)

protected = Unprotect[NonPositive];

NonPositive[ x_ /; Sign[x]<0 ] := True;
NonPositive[ x_ /; Sign[x]>0 ] := False;
NonPositive[ x_?Positive ] := False;
NonPositive[ x_ /; !Positive[x] ] := True;

Protect[protected];

(******************* Negative ********************************)

Unprotect[Negative];

Negative[ x_ /; Sign[x]<0 ] := True;
Negative[ x_ /; Sign[x]>0 ] := False;
Negative[ x_?NonNegative ] := False;
Negative[ x_ /; !NonNegative[x] ] := True;

Protect[Negative];

(******************* IntegerQ, EvenQ, OddQ ***************************)

Unprotect[IntegerQ,EvenQ,OddQ];

IntegerQ[x_Plus /;  (And @@ (IntegerQ /@ (List @@ x)))] := True;
IntegerQ[x_Times /; (And @@ (IntegerQ /@ (List @@ x)))] := True;
IntegerQ[(a_?IntegerQ) ^ ((b_?IntegerQ)?NonNegative)] := True;

EvenQ[x_Plus /;  (And @@ ((EvenQ[#] || OddQ[#] &) /@ (List @@ x)))] := 
  Not[Xor @@ (OddQ /@ (List @@ x)) ];
EvenQ[x_Times /; (And @@ (IntegerQ /@ (List @@ x))) &&
                 MemberQ[ EvenQ /@ (List @@ x), True ]] := True;
EvenQ[x_Times /; (And @@ (OddQ /@ (List @@ x)))] := False;
EvenQ[(a_?EvenQ) ^ ((b_?IntegerQ)?Positive)]    := True;
EvenQ[(a_?OddQ)  ^ ((b_?IntegerQ)?NonNegative)] := False;

OddQ[x_Plus /;  (And @@ ((EvenQ[#] || OddQ[#] &) /@ (List @@ x)))] := 
  Xor @@ (OddQ /@ (List @@ x));
OddQ[x_Times /; (And @@ (IntegerQ /@ (List @@ x))) &&
                 MemberQ[ EvenQ /@ (List @@ x), True ]] := False;
OddQ[x_Times /; (And @@ (OddQ /@ (List @@ x)))] := True;
OddQ[(a_?OddQ)  ^ ((b_?IntegerQ)?NonNegative)] := True;
OddQ[(a_?EvenQ) ^ ((b_?IntegerQ)?Positive)]    := False;

Protect[IntegerQ,EvenQ,OddQ];

(******************* InsertIndices ********************************)

(* What to do with inserted indices?  CovD or TCompute? 
   If indices are inserted into an expression of rank 0, first
   call TCompute to make sure it's a Component expression; 
   then call CovD to compute covariant derivatives. *)

(* First pull apart sums and scalar multiples. *)

InsertIndices[x_Plus, {i___}] := InsertIndices[#,{i}]& /@ x;

InsertIndices[f_ x_, {i__}] := 
  InsertIndices[f,{}] * InsertIndices[x, {i}] /;
    rank[f]===0 && rank[x]>0;

InsertIndices[x_, {}] := 
  Which[
    Head[rank[x]] =!= Integer,
      Print["ERROR: Invalid tensor expression ",x];
      ERROR[InsertIndices[x,{}]],
    rank[x] > 0,      x,
    ConstantQ[x],     x,
    True,             TCompute[x,{}]
    ];

InsertIndices[exp_, {i__}] := Module[ {rk},
  rk = rank[exp];
  Which[
    Head[rk] =!= Integer,
      Print["ERROR: Invalid tensor expression ",exp];
      ERROR[InsertIndices[exp,{i}]],
    rk === 0,
      CovD[ InsertIndices[exp,{}], {i} ],
    rk < Length[{i}],
      Message[Index::error, exp];
      ERROR[InsertIndices[exp,{i}]],
    rk > Length[{i}] && 
       !MemberQ[{Tensor,LeviCivitaConnection,Curvature}, Head[exp]],
      Message[Index::error, exp];
      ERROR[InsertIndices[exp,{i}]],
    !IndexQ[i],
      Print["ERROR: Invalid indices: ",{i}];
      ERROR[InsertIndices[exp,{i}]],
    !bundleOK[exp,{i}],
      0,
    True,
      TCompute[exp, {i}]
    ]
  ];    

(* bundleOK[ expression, {new indices} ] returns True if the bundles
   associated with new indices are correct for the corresponding slots
   of expression.  Because the indices may be complex, we have to be
   careful about the variances.  E.g. if slot 1 is a Covariant slot
   associated with complex bundle E, then you can insert a lower E
   index or an upper Conjugate[E] index, but not vice-versa. *)

bundleOK[exp_,{j___}] := Module[ { indexbunlist,tensorbunlist,
                                   indexvarlist,tensorvarlist,
                                   lists},
  indexbunlist  = Bundle /@ {j};
  indexvarlist  = Variance /@ {j};
  tensorbunlist = Take[Bundles[exp],Length[{j}]];
  tensorvarlist = Take[Variance[exp],Length[{j}]];
  If[MemberQ[indexbunlist,Null], Message[Bundle::error, 
    {j} [[ Position[indexbunlist,Null] [[1,1]] ]] ] ];
  lists = Transpose[{tensorbunlist,indexbunlist,tensorvarlist,indexvarlist}] //. 
     { {{Any},_,_,_} :> {{Any},Any,Covariant,Covariant},
       {buns_,index_,Covariant,Contravariant} :> 
         {buns,Conjugate[index],Covariant,Covariant},
       {buns_,index_,Contravariant,Covariant} :> 
         {buns,Conjugate[index],Covariant,Covariant} };
  Return[And @@ ( (MemberQ @@ (#[[{1,2}]]) &) /@ lists )]
  ];

(********************* Implicit index insertion *********************)

(* Here's where all the functions that can take indices inserted
   implicitly in brackets are defined *)

(* Explicit numbers *)
L/: n_Integer [i_L,j___] := InsertIndices[n,{i,j}];
U/: n_Integer [i_U,j___] := InsertIndices[n,{i,j}];
L/: n_Rational[i_L,j___] := InsertIndices[n,{i,j}];
U/: n_Rational[i_U,j___] := InsertIndices[n,{i,j}];
L/: n_Real    [i_L,j___] := InsertIndices[n,{i,j}];
U/: n_Real    [i_U,j___] := InsertIndices[n,{i,j}];
L/: n_Complex [i_L,j___] := InsertIndices[n,{i,j}];
U/: n_Complex [i_U,j___] := InsertIndices[n,{i,j}];

(* Tensors *)
HoldPattern[Tensor[x__][i___?IndexQ]] := 
  InsertIndices[Tensor[x],{i}];

(* Other Ricci functions. *)
HoldPattern[Curvature[x_][i___]]      := InsertIndices[Curvature[x],{i}];
HoldPattern[Del[x__][i___]]           := InsertIndices[Del[x],{i}];
HoldPattern[Div[x__][i___]]           := InsertIndices[Div[x],{i}];
HoldPattern[Extd[x_][i___]]           := InsertIndices[Extd[x],{i}];
HoldPattern[ExtdStar[x__][i___]]      := InsertIndices[ExtdStar[x],{i}];
HoldPattern[Grad[x__][i___]]          := InsertIndices[Grad[x],{i}];
HoldPattern[HodgeInner[x__][i___]]    := InsertIndices[HodgeInner[x],{i}];
HoldPattern[Int[x__][i___]]           := InsertIndices[Int[x],{i}];
HoldPattern[LeviCivitaConnection[x__][i___]] := 
                           InsertIndices[LeviCivitaConnection[x],{i}];
HoldPattern[Lie[x__][i___]]           := InsertIndices[Lie[x],{i}];
HoldPattern[RiemannTensor[x__][i___]] := InsertIndices[RiemannTensor[x],{i}];
HoldPattern[RicciTensor[x__][i___]]   := InsertIndices[RicciTensor[x],{i}];
HoldPattern[ScalarCurv[x__][i___]]    := InsertIndices[ScalarCurv[x],{i}];
HoldPattern[TensorProduct[x__][i___]] := InsertIndices[TensorProduct[x],{i}];
HoldPattern[Summation[x__][i___]]     := InsertIndices[x,{i}];
HoldPattern[Wedge[x__][i___]]         := InsertIndices[Wedge[x],{i}];
HoldPattern[Sym[x__][i___]]           := InsertIndices[Sym[x],{i}];
HoldPattern[Alt[x__][i___]]           := InsertIndices[Alt[x],{i}];

(* System functions *)
Unprotect[{Tr,Times,Power,Plus,Inner,Dot,Inverse,Transpose,Det}];
HoldPattern[Tr[x__][i___]]      := InsertIndices[Tr[x],{i}] /; IndexQ[i];
HoldPattern[Times[x__][i___]]   := InsertIndices[Times[x],{i}] /; IndexQ[i];
HoldPattern[Power[x__][i___]]   := InsertIndices[Power[x],{i}] /; IndexQ[i];
HoldPattern[Plus[x__][i___]]    := InsertIndices[Plus[x],{i}]  /; IndexQ[i];
HoldPattern[Inner[x__][i___]]   := InsertIndices[Inner[x],{i}] /; IndexQ[i];
HoldPattern[Inverse[x__][i___]] := InsertIndices[Inverse[x],{i}] /; IndexQ[i];
HoldPattern[Dot[x__][i___]]     := InsertIndices[Dot[x],{i}]   /; IndexQ[i];
HoldPattern[Transpose[x__][i___]] := InsertIndices[Transpose[x],{i}] /; IndexQ[i];
HoldPattern[Det[x__][i___]]     := InsertIndices[Det[x],{i}]   /; IndexQ[i];
Protect[{Times,Power,Plus,Inner,Dot,Inverse,Transpose,Det}];

(* Call DefineMathFunction so that expressions such as f[g][L[i]] will
   be interpreted correctly, where f is a mathematical function *)

DefineMathFunction[ {Sin,Cos,Tan,Sec,Csc,Cot}, 
                    Type -> Automatic, Quiet->True ];
DefineMathFunction[ Log, Quiet->True, Type->Complex ];

(* Allow implicit covariant derivatives of all Mathematica constants *)

Unprotect[{Catalan, Degree, E, EulerGamma, GoldenRatio, Pi}];
Catalan/:     Catalan    [i:(_L|_U)...] := InsertIndices[Catalan,{i}];
Degree/:      Degree     [i:(_L|_U)...] := InsertIndices[Degree,{i}];
E/:           E          [i:(_L|_U)...] := InsertIndices[E,{i}];
EulerGamma/:  EulerGamma [i:(_L|_U)...] := InsertIndices[EulerGamma,{i}];
GoldenRatio/: GoldenRatio[i:(_L|_U)...] := InsertIndices[GoldenRatio,{i}];
Pi/:          Pi         [i:(_L|_U)...] := InsertIndices[Pi,{i}];
Protect[{Catalan, Degree, E, EulerGamma, GoldenRatio, Pi}];

(*************************** TensorSimplify ************************)

(***************** defineTermwise, defineSumwise **********************)

(* The internal subroutine defineTermwise sets up the definitions that
   allow the simplification functions to work on a single term or
   specified terms of an expression, and also allow them to be
   threaded through lists of expressions. *)

defineTermwise[func_] := (
  func[ expr_List,  n___ ] := 
    func[ #, n ]& /@ expr;
  func[ expr_,      n_Integer,   opts___Rule ] := 
    func[ expr, {n}, opts ];
  func[ expr_Plus, {n__Integer}, opts___Rule ] := 
    MapAt[ func[#,opts]&, expr, List /@ {n}];
  func[ expr_, {1}, opts___Rule ] := 
    func[ expr, opts ];
  func[ expr_, {n__Integer}, opts___Rule ] := (
    Message[ Ricci::oneterm, expr ];
    Return[expr] );
  func[ expr_, opts__Rule ] := expr /;
    !optionsOK[ func, opts ];
  func[ expr_, (x_ /; Head[x]=!=Rule), opts___ ] := (
    Message[ Ricci::termspec, x ];
    Return[expr] );
  func[ ] := (
    Message[ Ricci::argm, func, 0, 1];
    Return[Null] );
  func[ expr_Plus, opts___Rule ] := 
    (func[#,opts]&) /@ expr;
  );

(* defineSumwise is similar to defineTermwise, except that it is for
   functions like CollectConstants or TensorSimplify that have to work
   on whole sums or parts of sums.  *)

defineSumwise[func_] := (
  func[ expr_List,  n___ ] := 
    func[ #, n ]& /@ expr;
  func[ expr_,      n_Integer,   opts___Rule ] := 
    func[ expr, {n}, opts ];
  func[ expr_Plus, {n__Integer}, opts___Rule ] := 
    Module[{others,i},
      others = Complement[Table[i,{i,Length[expr]}],{n}];
      func[ expr[[{n}]], opts ] + expr[[others]]
    ];
  func[ expr_, {1}, opts___Rule ] := 
    func[ expr, opts ];
  func[ expr_, {n__Integer}, opts___Rule ] := (
    Message[ Ricci::oneterm, expr ];
    Return[expr] );
  func[ expr_, opts__Rule ] := expr /;
    !optionsOK[ func, opts ];
  func[ expr_, (x_ /; Head[x]=!=Rule), opts___ ] := (
    Message[ Ricci::termspec, x ];
    Return[expr] );
  func[ ] := (
    Message[ Ricci::argm, func, 0, 1];
    Return[expr] );
  );

(*************************** TensorSimplify **************************)

defineSumwise[TensorSimplify];

TensorSimplify[expr_] :=
  expr //
  (CorrectAllVariances[#,Mode->OneDims]&)//(* put indices at correct altitude*)
  TensorExpand //        (* expand products & powers *)
  (AbsorbMetrics[#,Mode->NoOneDims]&) // (* eat all but 1-dim metrics *)
  PowerSimplify //       (* simplify bases & exponents *)
  RenameDummy //         (* rename dummy indices *)
  OrderDummy //          (* rearrange dummy names *)
  AbsorbMetrics //       (* absorb the remaining 1-dim metrics *)
  CollectConstants;      (* collect terms with same tensor factor *)

(*************************** SuperSimplify **************************)

defineSumwise[SuperSimplify];

SuperSimplify[expr_] :=
  expr //
  (CorrectAllVariances[#,Mode->OneDims]&) // 
  TensorExpand //
  (AbsorbMetrics[#,Mode->NoOneDims]&) //
  PowerSimplify //
  RenameDummy //
  (OrderDummy[#,Method->2]&) //
  AbsorbMetrics // 
  CollectConstants; 

(*************************** CovDSimplify **************************)

defineSumwise[CovDSimplify];

CovDSimplify[expr_] :=
  expr //
  TensorSimplify //
  OrderCovD //
  TensorSimplify;

(*************************** CollectConstants **************************)

defineSumwise[CollectConstants];

CollectConstants[expr_Plus] := Module[ {listofpairs},
  listofpairs = Sort[ getfactors /@ (List @@ expr)] //.
    {a___,{b_,c1_},{b_,c2_},d___} :> {a,{b,c1+c2},d};
  listofpairs = MapAt[simpconst,#,{2}]& /@ listofpairs;
             (*simplify constants while we have them!*)
  Return[Plus @@ ((Times @@ #)& /@ listofpairs)]
  ];

(* if not a sum, just simplify the constant factor. *)

CollectConstants[expr_] :=
  Times @@ MapAt[simpconst, getfactors[expr], {2}];

(* The following rule factors out a common I in the constant factor.
   Otherwise the presence of I will cause Mathematica to factor into
   Gaussian integers.
   We want (I m^2 + I) -> I(m^2 + 1), not (1 + I m)(I + m). *)

HoldPattern[ simpconst[ a:Plus[ (Times[Complex[0,_],_.]).. ] ] ] :=
  I simpconst[Expand[-I a]];

simpconst[x_] := If[LeafCount[x] > 50,
                    (*then*) Together[x],
                    (*else*) Factor[Together[x]]];

getfactors[a_Times] := getfactors /@ a;
getfactors[a_] := {1,a} /; ConstantQ[a];
getfactors[a_] := {a,1};
TensorFactor[a_] := First @ getfactors[a];
ConstantFactor[a_] := Last @ getfactors[a];

(************** FactorConstants ************************************)

defineTermwise[FactorConstants];

FactorConstants[x_] :=
  Times @@ MapAt[Factor, getfactors[x], {2}];

(************** SimplifyConstants **********************************)

defineTermwise[SimplifyConstants];

SimplifyConstants[x_] :=
  Times @@ MapAt[Simplify, getfactors[x], {2}];

(*************************** TensorExpand **************************)

(* This function does the same job as Expand on a tensor
   expression, but does not expand constant factors.  It also
   renames dummy indices as necessary to prevent conflicts.  *)

defineTermwise[TensorExpand];

TensorExpand[exp_] := Module[{newexp},
  newexp = exp //. a_Plus?ConstantQ :> fakeplus @@ a;
  newexp = newexp //. {
     HoldPattern[ (a_Times /; MemberQ[a,_Plus]) ] :> 
       Distribute[ a, Plus ],
     HoldPattern[ (a_Plus?(!ConstantQ[#]&))^(p_Integer?Positive) ] :>
       (a^(p-1) * # &) /@ NewDummy[a],
     HoldPattern[ (a_ /; GetDummyIndices[a] =!= {}) 
                ^ (p_Integer?Positive) ] :> 
       a^(p-1) * NewDummy[ a^1 ] 
     };
  Return[ newexp //. fakeplus -> Plus ]
  ];

(*************************** PowerSimplify **************************)

defineTermwise[PowerSimplify];

PowerSimplify[expr_] := expr //. 
  {head_[ x___, HoldPattern[Summation[y_]], z___ ] :> 
     head[x,y,z] /; head=!=Power,
   (a_ /; !FreeQ[a,L|U]) ^ p_Integer * (b_ /; !FreeQ[b,L|U]) ^ p_Integer :>
     (a b)^p,
   HoldPattern[ (a_ ^ b_) ^ c_ ] :> a ^ (b c)  /;
     IntegerQ[c] || (NonNegative[a] && Conjugate[b]===b),
   HoldPattern[ (a_ ^ b_Integer?Negative) ^ c_ ] :> (a^(-b))^(-c) } //.
  {HoldPattern[ (a_ ^ b_Integer) ^ c_ * (a_ ^ d_Integer) ] :> 
     (a^b) ^ (c + Quotient[d,b]) * (a ^ Mod[d,b]) /;
     Quotient[d,b] =!= 0} //.
  {Power[a_,b_] :>
    CollectConstants[TensorExpand[a]] ^ 
      CollectConstants[TensorExpand[b]] /;
      !IntegerQ[b] || b < 0
  };

(********************** CorrectAllVariances ****************************)

defineTermwise[CorrectAllVariances];

Options[CorrectAllVariances] := {Mode -> All};

CorrectAllVariances[exp_,opts___] := Module[{mode},
  mode = Mode /. {opts} /. Options[CorrectAllVariances];
  Switch[ mode,
    OneDims, Return[exp /. {HoldPattern[t_Tensor]:>correcttensorvariance[t,1]} //.
                             {x:(_Del|_Div|_Grad|_Extd|_ExtdStar|_Lie) :> 
                             CorrectAllVariances /@ x}],
    All,     Return[exp /. {HoldPattern[t_Tensor] :> correcttensorvariance[t,2]}],
    _,       Message[CorrectAllVariances::invalid, mode, 
                     "value for the Mode option"];
             Return[exp] ]
  ];

correcttensorvariance[ HoldPattern[Tensor[t_,ilist_,jlist_]], k_ ] :=
  Module[{newilist,newjlist,imetrics,jmetrics},
    If[ t===Basis || MetricQ[t] || {i,j}==={{},{}}, 
        (*then*) Return[Tensor[t,ilist,jlist]]];
    {newilist,imetrics} = correctindices[ ilist,
                                          Take[Variance[t],Length[ilist]],k];
    {newjlist,jmetrics} = correctindices[ jlist,
                                          Table[Covariant,{Length[jlist]}],k];
    Return[ Tensor[t,newilist,newjlist] * Times @@ Join[imetrics,jmetrics] ]; 
  ];

correctindices[ {},{},_ ] := {{},{}};

correctindices[ indexlist_, varlist_, k_ ] :=
  Transpose[ correctoneindex[#,k]& /@ Transpose[{indexlist,varlist}] ];

correctoneindex[ {L[i_],Covariant},    _ ] := {L[i],1};
correctoneindex[ {U[i_],Contravariant},_ ] := {U[i],1};

correctoneindex[ {L[i_],Contravariant}, k_ ] := 
  If[ k===1 && Dimension[Bundle[i]] =!= 1,
      (*then*) {L[i],1},
      (*else*) {U[#], InsertIndices[ Metric[Bundle[i]], {L[i],L[#]} ]}& @
                 NewBundleSymbol[Conjugate[Bundle[i]]] ];

correctoneindex[ {U[i_],Covariant}, k_ ] := 
  If[ k===1 && Dimension[Bundle[i]] =!= 1,
      (*then*) {U[i],1},
      (*else*) {L[#], InsertIndices[ Metric[Bundle[i]], {U[i],U[#]} ]}& @
                 NewBundleSymbol[Conjugate[Bundle[i]]] ];

(********************** Expand **************************************)

Unprotect[Expand];
Expand[x_] := TensorExpand[x] /; !FreeQ[x, (L|U)];
Protect[Expand];

(*************************** NewDummy **************************)

(* This function replaces all dummy index pairs in an expresion by
   unique dummy names *)

defineTermwise[NewDummy];

NewDummy[expr_] := Module[{allindices,pairs,newdums,rules},
  (* allindices = list of all the indices in expr *)
  allindices = GetIndices[expr];
  (** pairs = a list of the lower member of paired indices in expr.  **)
  pairs = Intersection[
              (Pair /@ Cases[allindices,U[_]]), Cases[allindices,L[_]]
          ] /. x_[Bar] :> x;
  (* newdums = list of unique names of the form dN *)
  newdums = NewBundleSymbol /@ Bundle /@ pairs;
  pairs = (pairs /. {L[i_] :> i});
  (* now create the rule which will replace paired indices with
     dummy pairs.
     step 1: pairs   = { i1, i2, ..., iN }
             newdums = { d(n+1),  d(n+2),     , d(n+N) }
             the "Transpose" command creates a list that looks like
             { {i1,d(n+1)}, ..., {iN,d(n+N)} }
     step 2: the dummyrule function transforms a list like {i,d}
             into a list of rules like
               {i -> d}
     step 3: the "Flatten" flattens out the list of rules.
  *)
  newdums = Select[
              Transpose[ {pairs, newdums} ],
              Last[#] =!= Null &
              ];
  rules = Flatten[ dummyrule /@ newdums ];
  (* Finally, apply the rules. *)
  Return[(expr /. rules)]
  ];

dummyrule[{i_,d_}] := {i -> d};

(*************************** RenameDummy *************************)

defineTermwise[RenameDummy];

RenameDummy[expr_] := Module[{oldfrees, olddums, tensorlist, newdums,
                              oldnewpairs},

   (* oldfrees,olddums = lists of free and dummy index names in expr. *)

   oldfrees = Union[ First /@ GetFreeIndices[Hold[expr]] /. {x_[Bar]->x} ];
   olddums  = Union[ First /@ GetDummyIndices[Hold[expr]] /. {x_[Bar]->x} ];

   (* If there are no dummy indices, there's nothing to do. *)

   If[ olddums==={}, Return[ expr ] ];

   (* get a list of the tensors associated with each index name (with
   indices blanked out), so we can canonicalize the initial choice
   of dummy names as much as possible.  Then sort olddums according 
   to their associated tensor names and positions. *)

   tensorlist = gettensors[ tenfactors, # ]& /@ olddums;
   {tensorlist,olddums} = 
     Transpose[ Sort[ Transpose[{tensorlist,olddums}] ] ];

   (* Change names of dummy indices *)

   Resetcanbundum[oldfrees];
   newdums = canbundum /@ olddums;

   (* oldnewpairs is a list like {{i1,d1},{i2,d2},...}, where i1=old
   dummy, d1=new dummy, etc.  (Ignore dummies with no bundle if they
   happen to sneak in somehow.)  *)

   oldnewpairs = Select[ Transpose[ {olddums, newdums} ],
                         Last[#] =!= Null & ];

   (* create the list of rules i1->d1 etc. & apply the rules *)

   Return[ expr /. Flatten[(Rule @@ # )& /@ oldnewpairs] ]
   ];

(*************************** OrderDummy **************************)

(* This is the heart of TensorSimplify.  It attempts to put the
   indices in an expression into a canonical form by rearranging the
   names of dummy indices.  *)

defineTermwise[OrderDummy];
Options[OrderDummy] :=
  {Method -> 1};

OrderDummy[expr_, opts___] := Module[{method, dums, tenfactors,
                              confactors, tensorlist, pairlist,
                              possiblechanges, rules, i, j},

   method = Method /. {opts} /. Options[OrderDummy];
   If[ !MemberQ[ {0,1,2}, method ], 
       (*then*) Message[OrderDummy::invalid, method, "Method option"];
                Return[expr] ];

   (* dums = lists of dummy index names in expr. *)

   dums  = Union[ First /@ GetDummyIndices[Hold[expr]] /. {x_[Bar]->x} ];

   (* If there are no dummy indices, just return. *)

   If[ dums==={}, Return[ expr ] ];

   (* First separate the constant factors from the tensor factors. *)

   {tenfactors,confactors} = getfactors[ expr ];

   (* Hold the tensor part of the expression, so we don't get into
   infinite loops with symmetry rules working automatically *)

   tenfactors = Hold[Evaluate[tenfactors]];

   (* get a list of the tensors associated with each index name. *)

   tensorlist = gettensors[ tenfactors, # ]& /@ dums;

   (* Now correct "slopes" of dummy index pairs by interchanging upper
   and lower members of pairs.  "pairlist" is a list of pairs of
   positions where dummy pairs occur.  *)

   pairlist = Flatten[ getpairlist[tenfactors,#]& /@ dums, 1 ];

   {tenfactors,pairlist} = correctvarianceandbars[tenfactors,pairlist];

   (* Next we try interchanging each of the remaining pairs to see if
   the expression looks better that way.  When we're done, we need no longer
   hold tenfactors. *)

   tenfactors = ReleaseHold[fixslopes[tenfactors,pairlist]];

   (* Finally, we're going to try permuting the names of the dummy indices. 
      The "Method" option controls how hard we work. *)

   Switch[ method,

     0, (* don't change names at all *)

       Null,

     1, (* normal level of simplification: only try switching index names
        a pair at a time. *)

        (* Make a list of pairs of indices whose names we might want to
        interchange.  They must have the same bundle, and at least one
        tensor in common. *)

        possiblechanges = 
          Select[ Flatten[ Table[{i,j}, {i,Length[dums]-1},
                                        {j,i+1,Length[dums]}], 
                           1],
                  Bundle[dums[[ #[[1]] ]]]===Bundle[dums[[ #[[2]] ]]] &&
                  Intersection[ tensorlist[[ #[[1]],2 ]], 
                                tensorlist[[ #[[2]],2 ]] ] =!= {} &
                ];

        rules = (exchangerule[ dums[[#]] ] &) /@ possiblechanges;

        (* Now try changing the names corresponding to each pair.  Keep
        the better expression.  From this point on, we don't need to
        hold tenfactors any more. *)

        tenfactors = fixnames[ tenfactors, rules ],

     2, (* work harder to simplify the expression by trying
        all allowable permutations of the index names *)

        possiblechanges = (List /@ Range[Length[dums]]) //. 
                          {a___, {b__}, {c_}, d___} :> {a, {b,c}, d} /; 
                          SameQ @@ (Bundle /@ dums[[{b,c}]]);

        possiblechanges = shuffles[possiblechanges];

        rules = permuterule[ dums, dums[[#]] ]& /@ possiblechanges;

        tenfactors = best[ (tenfactors /. # &) /@ rules ]

     ];

   (* Finally, construct the expression to return to the caller. *)  

   Return[ confactors * ReleaseHold[tenfactors] ];

   ];

(******************** Subroutines for OrderDummy ********************)

(* SortOn[listoflists,n] sorts "listoflists" by its nth element.
   SortOn[listoflists,{n1,n2,...}] sorts by elements n1,n2,... in order. *)

SortOn[list_,n_] := Last /@ Sort[ {#[[n]],#}& /@ list ];

(* gettensors[exp,i] compiles a list of tensors and positions in which
   index i occurs in exp (with all indices replaced by Null).  This
   enables us to order our initial choice of dummy index names
   according to the tensor "patterns" in which indices occur, which is
   independent of the existing names. *)

gettensors[exp_,i_] := Module[{ipos},
  ipos = Position[ exp, L[i] | L[i[Bar]] | U[i] | U[i[Bar]] ];
  Return[ 
    {Bundle[i],
     Sort[(HeldPartAt[exp,Drop[#,-2]]& /@ ipos) //. 
          {L[_] -> Null, U[_] -> Null} ],
     ipos}
    ];
  ];

(* getpairlist[exp,i] returns a list of pairs of positions
   showing where pairs of dummy indices occur.  If a given occurrence of L[i]
   is paired with more than one U[i], ignore all but the first.  Multiple
   occurrences should happen only with one-dimensional bundles, which we're
   not considering here. *)

getpairlist[exp_,i_] := Module[{productops,lowpos,highpos,pairpos},
  productops = {Times,TensorProduct,Wedge,Int,Inner,HodgeInner,Dot};
  lowpos = Position[ exp, L[i] | L[i[Bar]] ];
  highpos = pairloc[exp,#,productops]& /@ lowpos;
  (* make a list of pairs, throwing out occurrences of L[i] with no U[i] *)
  pairpos = Cases[ Transpose[{lowpos,highpos}], {_,{_,___}} ];
  Return[ { #[[1]], First[ #[[2]] ] }& /@ pairpos ]
  ];

(* pairloc[exp, position, validheads] picks out the locations in exp
   of the paired dummy index for the index located at "position".  We
   only look at stuff combined the operator heads listed in
   "validheads".  The function value is a list of positions in "exp"
   where valid occurrences of the paired index lie. *)

pairloc[ exp_, pos_, validheads_ ] := 
  Module[ {heads,goods,subexppos,subexp,pairs},
  (* get the chain of heads leading to the index at position "pos" 
     but drop the final Tensor, List, & L or U *)
  heads = Drop[ headlist[ exp, pos ], -3 ];
  (* figure out which heads are legitimate pair combiners *)
  goods = MemberQ[ validheads, # ] & /@ heads;
  (* get the subexpression in which to look *)
  subexppos = Take[ pos, Last[ Position[ goods, False ] ][[1]] ];
  subexp = HeldPartAt[exp,subexppos];
  pairs = Position[ subexp, Pair[ PartAt[exp,pos] ] ];
  (* if there are none, give up. *)
  If[ pairs==={}, Return[{}] ];
  (* now throw out those that are embedded in non-product operators *)
  heads = Take[ headlist[subexp,#], {2,-4} ]& /@ pairs;
  goods = Map[ MemberQ[ validheads, # ]&, heads, {2} ];
  pairs = Cases[ Transpose[{pairs,goods}], 
                 {{___},{True...}} | {{___},{}} ];
  If[ pairs =!= {},
      pairs = First @ Transpose @ pairs ];
  Return[ Join[subexppos,Drop[#,1]]& /@ pairs ]
  ];

(* headlist[expression,position] gives the chain of heads leading to 
   the object at "position" *)

headlist[ exp_, pos_ ] := Module[{i},
  PartAt[ exp, # ]& /@
    Table[ Append[ Take[pos,i], 0 ], {i,0,Length[pos]} ]
  ];

(* correctvarianceandbars[exp,pairlist] looks at each of the dummy pairs
   in pairlist, and first tries to put both at "natural" altitude.  If
   that's not possible, it then unbars any complex pairs.  When finished,
   any pairs that are at natural altitude and all complex pairs are removed
   from the list.  Indices that appear with Basis are removed from the
   list untouched. *)

correctvarianceandbars[ exp_, pairlist_ ] :=
  Module[{newexp,newpairlist,j,pair,ind,bun,ten,variance,bar},
  newexp = exp; 
  newpairlist = pairlist;
  Do[ (* j = 1 to Length[newpairlist] *)
    pair = newpairlist[[j]];
    {ind,bun,ten,variance,bar} = 
      Transpose[ getindexdata[exp,#]& /@ pair ];
    Which[
      ten[[1]] === Basis || ten[[2]] === Basis,
        (* If either tensor appears with Basis, don't change it. *)
            newpairlist[[j]] = Null,
      variance[[{1,2}]] === {False,False},
        (* Incorrect variance: this pair should be swapped *)
            newexp = doswap[newexp,pair];
            (* Don't swap this pair again *)
            newpairlist[[j]] = Null,
      variance[[{1,2}]] === {True,True},
        (* Already correct variance: this pair should be left in peace *)
            newpairlist[[j]] = Null,
      bar[[1]]===Bar,
        (* Barred pair that can be unbarred *)
            newexp = doswap[newexp,pair];
            newpairlist[[j]] = Null,
      Type[ bun[[1]] ] === Complex,
        (* Unbarred complex pair that shouldn't be swapped *)
            newpairlist[[j]] = Null ],
    {j,Length[newpairlist]} ];
  Return[{newexp, Select[ newpairlist, # =!= Null & ]}] 
  ];

(* getindexdata returns an "indexdata" list about a single index position:
   1. the position
   2. the bundle
   3. the tensor name
   4. True if correct variance, False if not
   5. Bar if barred, Null if not  *)

getindexdata[ exp_, pos_ ] := Module[{ten,ind,indexpos},
  ten = HeldPartAt[ exp, Drop[pos,-2] ] [[1,1]];  (* the tensor name *)
  ind = PartAt[ exp, pos ];           (* the index *)
  (* indexpos = the slot in ten in which ind occurs *)
  indexpos = If[ pos[[-2]]===3,
       pos[[-1]]+Rank[ten],  (* differentiated index *)
       pos[[-1]] ];          (* ordinary index *)
  Return[
    { pos,                (* the position *)
      Bundle[ind],        (* the bundle *)
      ten,                (* the tensor name *)
      var[ten,indexpos]===var[ind],  (* True if correct variance *)
      If[ FreeQ[ind,Bar], 
          Null,  (* unbarred index *)
          Bar ]  (* barred index *)
    } ]
  ];

(* var[tensor,pos] gets the natural variance of the index at
   position pos in tensor; var[ind] returns the actual variance
   of an index. *)

var[t_[Bar], n_] := var[t,n];

var[t_, n_] :=
  If[ n > rank[t],
    (*then*) Covariant,       (* Differentiated indices are Covariant *)
    (*else*) Variance[t][[n]] (* The others are what they are. *)
   ];

var[L[_]] := Covariant;
var[U[_]] := Contravariant;

(* fixslopes makes one pass through the list of pairs to
   see if the expression looks better with opposite slope *)

fixslopes[ exp_, pairlist_ ] := 
  Return[ First @ Fold[ pickbestslope, 
                        {exp,fixforcomparison[exp]}, 
                        pairlist ] ];

(* pickbestslope[ {x, fixedx}, {pos1,pos2} ] tries changing the slope
   of the indices at locations pos1, pos2, to see if the expression
   improves (in lexical order).  We assume fixedexpression is the
   result of sending expression through fixforcomparison.  The
   returned value is either {x,fixedx} or the revised version. *)

pickbestslope[ {x_,fixedx_}, pair_ ] := Module[{newx,newfixedx},
  If[ x===Hold[0], Return[{Hold[0],Hold[0]}] ];
  newx = doswap[x,pair];
  newfixedx = fixforcomparison[newx];
  Switch[ better[fixedx,newfixedx],
          True,  Return[{x,fixedx}],
          False, Return[{newx,newfixedx}],
          0,     Return[{Hold[0],Hold[0]}]
        ];
  ];

fixnames[ exp_, rules_ ] := 
  Return[ 
    First @ Fold[ pickbestnames, 
                  {exp, fixforcomparison[exp]}, 
                  rules ] ];

pickbestnames[ {x_, fixedx_}, rule_ ] := Module[{newx,newfixedx},
  If[ x===0, Return[{0,Hold[0]}] ];
  newx = x /. rule;
  newfixedx = fixforcomparison[newx];
  Switch[ better[fixedx,newfixedx],
          True,  Return[{x,fixedx}],
          False, Return[{newx,newfixedx}],
          0,     Return[{0,Hold[0]}]
        ];
  ];

(* doswap actually raises and lowers a specific pair of indices.  Note
   that the conjugate of a Real index is itself. *)

doswap[exp_,pair_] := Module[{newexp,ind1,ind2},
  ind1 = PartAt[ exp, pair[[1]] ];
  ind2 = PartAt[ exp, pair[[2]] ];
  newexp = ReplacePart[ exp,
                        Pair[Conjugate[ ind1 ]],
                        pair[[1]] ];
  newexp = ReplacePart[ newexp,
                        Pair[Conjugate[ ind2 ]],
                        pair[[2]] ];
  Return[newexp]
  ];

(* fixforcomparison[x] reevaluates x once, then replaces Times, Plus,
   and Power by dummy operations, because Mathematica has weird rules
   for ordering polynomials.  It also gloms all differentiated indices
   together with the undifferentiated ones, so that in effect indices
   are ordered as if the semicolon weren't there.  Finally, any
   overall minus signs are pulled out in front of the Hold, so they
   can be ignored when lexically comparing later. *)

fixforcomparison[ exp_ ] := 
    Hold[Evaluate[ReleaseHold[exp]]] //. 
                 {Times -> times, Power -> power, Plus -> plus,
                  Tensor[t_,{i___},{j___}] :> Tensor[t,{i,j}]} //.
                 {times[x___,-1,y___] :> - times[x,y],
                  times[x_] :> x};

(* better[ exp1, exp2 ] returns True if exp1 is lexically smaller than
   exp2 (after stripping signs), False if the other way around, and 0
   if they differ exactly by a sign. *)

better[ Hold[- x_], Hold[x_] ] := 0;
better[ Hold[x_], Hold[- x_] ] := 0;
better[ Hold[- x_], Hold[- y_] ] := OrderedQ[{Hold[x],Hold[y]}];
better[ Hold[- x_], Hold[y_] ]   := OrderedQ[{Hold[x],Hold[y]}];
better[ Hold[x_], Hold[- y_] ]   := OrderedQ[{Hold[x],Hold[y]}];
better[ Hold[x_], Hold[y_] ]     := OrderedQ[{Hold[x],Hold[y]}];

(* shuffles[ { {a,b}, {c,d,e} } ] returns a list of all permutations of
   {a,b,c,d,e} which consist of permutations of {a,b} and {c,d,e} 
   separately.  a,b,c,d,e must be integers.  *)

shuffles[{}] := {};
shuffles[ listoflists_ ] := Module[{list},
  list = Permutations /@ listoflists /. {i__Integer} :> sublist[i];
  Return[ Flatten[ Outer @@ Prepend[ list, Join ] ] /. sublist -> List ]
  ];

(* best[ {exp1, exp2, ... } ] returns the expression that is lexically
   smallest.  Signs are stripped off when comparing; if any two
   expressions differ by a sign, then best returns zero.  We assume
   that fixforcomparison has already been applied. *)

best[ explist__ ] := Module[{newlist},
  newlist = fixforcomparison /@ explist;
  newlist = (Replace[ #, -x_ :> x ]&) /@ newlist;
  newlist = Sort[ Transpose[ {newlist,explist} ] ];
  If[ MatchQ[ newlist,
              {___,{x_,-y_},{x_,y_},___}  ],
      (*then*) Return[0] ];
  Return[ newlist[[1,2]] ];
  ];

(* exchangerule takes a pair of objects and returns a set of rules that
   interchanges them. *)

exchangerule[{a_,b_}] := {a -> b, b -> a};

(* permuterule takes two lists, and returns a set of rules that transforms
   each item in the first list to the corresponding item in the second *)

permuterule[x_List,y_List] := (Rule @@ #)& /@ Transpose[{x,y}];

(************************* canbundum ************************************)

(* canbundum[index] returns the next dummy index name for index's
   bundle in canonical order.  The order is: first go through the
   user's list of indices associated with the bundle in the order
   given; when that list is used up, start generating new indices of
   the form i1, i2, etc., using the last index name in the list as
   prefix.  However, make sure as we go that we do not duplicate any
   indices in "freeindexlist".

   Before calling canbundum, we have to call Resetcanbundum with a
   list of free index names to avoid.  Because the list of free
   indices and index counters for each bundle are stored in global
   variables, OrderDummy is NOT recursive.  *)

Resetcanbundum[freeindices_] := (
  Clear[canbundumcount];
  canbundumcount[_] = 0;
  Clear[indicesusedup];
  indicesusedup[_]  = False;
  freeindexlist = freeindices;
  );

canbundum[i_] := Module[{indexlist,bundle},
  bundle = Bundle[i];
  If[bundle===Null, Return[Null] ];
  canbundumcount[bundle] = canbundumcount[bundle] + 1;
  indexlist = Complement[BundleIndices[bundle],freeindexlist];
  If[!indicesusedup[bundle] &&
     canbundumcount[bundle] === Length[indexlist] + 1,
    (*then*) indicesusedup[bundle] = True;
             canbundumcount[bundle] = 1
    ];
  If[indicesusedup[bundle],
    (*then*)
      sym = cansymbol[bundle];
      If[!(Dimension[bundle] === 1),
        While[ MemberQ[freeindexlist,sym],
               canbundumcount[bundle] += 1;
               sym = cansymbol[bundle]
             ],sym];
      Return[sym],
    (*else*)
      Return[ indexlist [[ canbundumcount[bundle] ]] ]
    ]
  ];

cansymbol[bundle_] := Module[{i},
  i = ToExpression[
    StringJoin[
      ToString[BundleDummy[bundle]],
      If[Dimension[bundle]===1,
        "",
        ToString[ canbundumcount[bundle] ]
        ]
      ]
    ];
  newsymsub[bundle,i];   (* See Index.m for definition of "newsymsub". *)
  Return[i]
  ];

(*************************** TensorCancel ********************************)

defineTermwise[TensorCancel];

TensorCancel[exp_] := Module[{scalarfactors,tensorfactors},
  {scalarfactors,tensorfactors} = 
    getScalarFactors @ PowerSimplify @ AbsorbMetrics @ exp;
  scalarfactors = CollectConstants /@ OrderDummy /@ RenameDummy /@ 
                  scalarfactors;
  scalarfactors = scalarfactors //. {
    {a___,x_^p_.,b___,x_^q_.,d___} :> {a,x^(p+q),b,d},
    {a___,HoldPattern[Summation[x_]]^p_,b___,x_^q_.,d___} :> {a,x^(p+q),b,d},
    {a___,x_^q_.,b___,HoldPattern[Summation[x_]]^p_,d___} :> {a,x^(p+q),b,d}
    };
  Return[ OrderDummy @ RenameDummy[ 
          (Times @@ (NewDummy /@ scalarfactors)) * (Times @@ tensorfactors) ]
        ];
  ];

getScalarFactors [exp_Times] := Module[{factors,scalarfactors,tensorfactors},
  factors = (List @@ exp)//.
     {{a___,x_,b___,y_,c___} :> {a, b, c, x y} /;
     !(Intersection[GetIndices[x], Pair /@ GetIndices[y]] === {}) &&
     !FreeQ[x,Tensor[_,_,_]]
     };
  scalarfactors = Select[factors, GetFreeIndices[#]==={}&];
  tensorfactors = Complement[factors, scalarfactors];
  Return[{scalarfactors,tensorfactors}];
  ];

getScalarFactors[exp_] := {{exp},{}} /; GetFreeIndices[exp] === {};
getScalarFactors[exp_] := {{},{exp}};

(********************** AbsorbMetrics *******************************)

defineTermwise[AbsorbMetrics];

Options[AbsorbMetrics] := {Mode -> All};

AbsorbMetrics[ exp_, opts___ ] := Module[{heldexp, metrics, mode},
  mode = Mode /. {opts} /. Options[AbsorbMetrics];
  If[ !MemberQ[ {All,NoOneDims}, mode ], 
      Message[AbsorbMetrics::invalid, mode, "value for the Mode option"];
      Return[exp] ];
  heldexp = expandpowers[Hold[exp]];
  metrics = Position[ heldexp, _Tensor?MetricQ ];
  Return[ ReleaseHold[ Fold[ eatmetric[##,mode]&, heldexp, metrics ] ] ]
  ];

eatmetric[ exp_, pos_, mode_ ] := 
  Module[{metric,hds,poslist,newexp,newindex},
  metric = HeldPartAt[ exp, pos ];
  indlist = metric[[1,2]];
  (* If the metric doesn't have exactly two indices, or if the indices
     are already paired, ignore it. *)
  If[ Length[indlist] =!= 2 || PairQ @@ indlist,
      (*then*) Return[exp] 
    ];
  (* see if the metric is equivalent to Kronecker[L[i],U[j]].  This is
     the case whenever the two indices have different altitudes.
     (Since the expression is held, the metric may not have been
     replaced by Kronecker yet.)  If so, it can be paired even with
     indices inside differential ops. *)
  hds = If[ SameQ @@ Head /@ indlist,  
            (*then*) {Times,TensorProduct,Wedge,Int,Inner,HodgeInner,Dot},
            (*else*) {Times,TensorProduct,Wedge,Int,Inner,HodgeInner,Dot,
                      Del,Grad,Div,Extd,ExtdStar,Lie} ];
  (* get the lists of indices paired with the left or right index of
     our metric.  *) 
  poslist = Join[ pairloc[ exp, Join[pos,{2,2}], hds ], 
                  pairloc[ exp, Join[pos,{2,1}], hds ] ];
  (* Throw away any paired indices that appear with Basis 
     (unless the "metric" is Kronecker). *)
  If[ SameQ @@ Head /@ indlist,
      poslist = 
        Select[ poslist, 
                tensornameat[ exp, # ] =!= Basis & ] ];
  (* If the Mode option is NoOneDims (set by TensorSimplify), also throw
     away one-dimensional paired indices that don't appear with metrics. *)
  If[ mode === NoOneDims,
      poslist =
        Select[ poslist,
                Dimension[Bundle[PartAt[exp,#]]] =!= 1 ||
                MetricQ[ tensornameat[exp,#] ]& ] ];
  (* Now sort the list by (a) whether the variance is correct (wrong first)
                          (b) whether the index is barred (barred first)
                          (c) position of the index *)
  poslist = SortOn[ getindexdata[ exp, # ]& /@ poslist,
                    {4,5,1}];
  (* Finally, eat the index that comes up first in the sort, if there is one.
     Replace the metric by 1, to preserve the overall structure of the 
     expression, so the positions of other metrics will still be valid. *)
  If[ poslist === {},
      (*then*) newexp = exp,
      (*else*) newindex = 
                 If[ PartAt[ exp, poslist[[1,1]] ] === Pair[metric[[1,2,1]] ],
                     (*then paired with first index*) metric[[1,2,2]],
                     (*else paired with second *)     metric[[1,2,1]] ];
               newexp = ReplacePart[ exp, 
                                     newindex, 
                                     poslist[[1,1]] ];
               newexp = ReplacePart[ newexp, 1, pos ];
    ];
  Return[newexp];               
  ];

(* The following subroutine is used by AbsorbMetrics to replace powers
   of indexed tensor expressions by multiplication.  This is because
   metrics of 1-dimensional bundles may show up to a positive power. *)

expandpowers[ exp_ ] := Module[{powerloc},
  powerloc = Position[ exp, Power[ (x_ /; !FreeQ[x,_L|_U]),
                                   (y_Integer /; y > 0)] ];
  Return[ Fold[ expandonepower, exp, powerloc ] ];
  ];

expandonepower[ exp_, {pos___} ] := Module[{base,exponent},
  base = exp[[pos,1]];
  exponent = exp[[pos,2]];
  Return[ insertmultiplication[ exp, {pos}, Table[base,{exponent}] ] ];
  ];

insertmultiplication[ exp_, {pos___}, {factors__} ] := 
  ReplaceHeldPart[ exp, Hold[Times[factors]], {pos} ];

tensornameat[heldexp_, pos_] := 
  PartAt[heldexp, Append[ Drop[pos,-2], 1]];

(*********** LowerAllIndices, LowerIndicesRule **********************)

defineTermwise[LowerAllIndices];

LowerAllIndices[expr_] :=
  expr //.
  { Tensor[t_,x___,{i___,U[j_],k___},y___]
      :> (Tensor[t,x,{i,L[#],k},y] *
          InsertIndices[ Metric[Bundle[j]], {U[j],U[#]}] &)
            [ NewBundleSymbol[ Conjugate[Bundle[j]] ] ] /;
         !MetricQ[t] && t=!=Basis,
    HoldPattern[Tensor[g_,{a___,L[b_],c___},{}] *
         Tensor[g_,{i___,U[b_],j___},{}]] :>
           Tensor[g,{a,c,i,j},{}] /; MetricQ[g]
  };

LowerIndicesRule[Tensor[name_,_,_]] :=
  { Tensor[name,x___,{i___,U[j_],k___},y___]
    :> (Tensor[name,x,{i,L[#],k},y] *
        InsertIndices[ Metric[Bundle[j]], {U[j],U[#]} ] &)
            [ NewBundleSymbol[ Conjugate[Bundle[j]] ] ] /;
       !MetricQ[name] && name =!= Basis
  };

LowerIndicesRule[{}] := {};
LowerIndicesRule[{i___}] := LowerIndicesRule[i];
LowerIndicesRule[i_,j__] := Flatten[LowerIndicesRule /@ {i,j}];
LowerIndicesRule[i_Symbol] := LowerIndicesRule[ Union[ {U[i],U[i[Bar]]} ] ];

LowerIndicesRule[U[j_]] :=
  { Tensor[name_,x___,{i___,U[j],k___},y___] :>
       (Tensor[name,x,{i,L[#],k},y] *
        InsertIndices[ Metric[Bundle[j]], {U[j],U[#]} ] &)
           [ NewBundleSymbol[ Conjugate[Bundle[j]] ] ] /;
       !MetricQ[name] && name =!= Basis
  };

RaiseIndicesRule[L[j_]] :=
  { Tensor[name_,x___,{i___,L[j],k___},y___] :>
       (Tensor[name,x,{i,U[#],k},y] *
        InsertIndices[ Metric[Bundle[j]], {L[j],L[#]} ] &)
           [ NewBundleSymbol[ Conjugate[Bundle[j]] ] ] /;
       !MetricQ[name] && name =!= Basis
  };

(*********************** BasisExpand ********************************)

defineTermwise[BasisExpand];

BasisExpand[exp_Times] := 
  TensorExpand[BasisExpand /@ exp] /; !Rank[exp]===0;
BasisExpand[exp^(p_Integer)] := 
  TensorExpand[BasisExpand[exp]^p] /; !Rank[exp]===0;
BasisExpand[exp_Wedge] := BasisExpand /@ exp;
BasisExpand[exp_TensorProduct] := BasisExpand /@ exp;
BasisExpand[exp_Alt] := BasisExpand /@ exp;
BasisExpand[exp_Sym] := BasisExpand /@ exp;
BasisExpand[exp_] := TensorExpand[InsertIndices[exp,{}]] /; Rank[exp]===0;

BasisExpand[exp_] := Module[{producttype,list,indices,n},
  producttype = Which[
                 SymmetricQ[exp],     Times,
                 SkewQ[exp],          Wedge,
                 True,                TensorProduct];
  list = Transpose[{Bundles[exp],Variance[exp]}];
  If [!FreeQ[list,Null],
     Message[Bundle::error,exp];
     Return[ERROR[exp]]
     ];
  indices = NewIndex /@ list;
  TensorExpand[ Plus @@ Flatten [
    Outer @@ Join[
        {InsertIndices[exp,{##}] * producttype @@
          Table[Tensor[Basis,{Pair[ {##}[[n]] ]},{}], {n,Length[{##}]}] &},
        indices
        ]
    ] /
    If[producttype === Wedge,
      (*then*) WedgeFactor[Length[indices]],
      (*else*) 1]]
  ];

(************************** BasisGather *************************)

SetAttributes[makerule,HoldRest];

BasisGather[exp_, {tensors__}] :=
  TensorExpand[exp] //. 
    Flatten[{BasisGatherRule /@ {tensors},
             BasisGatherRule /@ (Times[-1,#]&) /@ {tensors}}];

BasisGather[exp_, tensor_] := BasisGather[exp, {tensor}];

BasisGatherRule[tensor_ /; Rank[tensor]===1 ] := Module[{target},
  target = targetform[tensor];
  Return[{
    makerule[ target, x * tensor ],
    makerule[ 
      FixedPoint[ 
        Replace[
          #, 
          {HoldPattern[_[a__ * Tensor[Basis,{i_},{}]]] :>
             HoldPattern[a*Wedge[y___,Tensor[Basis,{i},{}],z___]],
           HoldPattern[_[b__ + a__ * Tensor[Basis,{i_},{}]]] :>
             HoldPattern[b+a*Wedge[y___,Tensor[Basis,{i},{}],z___]]}]&,
        target],
      x * Wedge[y,tensor,z]],
    makerule[ 
      FixedPoint[ 
        Replace[
          #, 
          {HoldPattern[_[a__ * Tensor[Basis,{i_},{}]]] :>
             HoldPattern[a*TensorProduct[y___,Tensor[Basis,{i},{}],z___]],
           HoldPattern[_[b__ + a__ * Tensor[Basis,{i_},{}]]] :>
             HoldPattern[b+a*TensorProduct[y___,Tensor[Basis,{i},{}],z___]]}]&,
        target],
      x * TensorProduct[y,tensor,z]]
    }]];

BasisGatherRule[tensor_?SymmetricQ ] :=
  {makerule[ targetform[tensor], x * tensor ]};

BasisGatherRule[tensor_?SkewQ ] :=
  {makerule[ targetform[tensor] /.
                HoldPattern[Wedge[a___]] :>
                Wedge[y___,a,z___],
             x * Wedge[y,tensor,z]]};

BasisGatherRule[tensor_ ] :=
  {makerule[ targetform[tensor] /.
                HoldPattern[TensorProduct[a___]] :>
                TensorProduct[y___,a,z___],
             x * TensorProduct[y,tensor,z]]};

targetform[tensor_] := 
  HoldPattern[Evaluate[NewDummy[TensorExpand[x_. * BasisExpand[tensor]]]]];

makerule[ lhs_, rhs_ ] :=
  Module[ {indices,buns,bunpairs,patternrules,timesrules,w},
    indices = First /@ GetDummyIndices[lhs] /. v_[Bar]:>v;
    buns = Bundle /@ indices;
    bunpairs = Transpose[{indices,buns}];
    patternrules = (# -> Pattern[#,Blank[]]&) /@ indices;
    timesrules = (HoldPattern[Times[ a___,Times[b__],c___ ]] :> Times[a,b,c]);
    Return[ (lhs :> rhs /; And @@ (Bundle[First[#]]===Last[#]&) /@ w)
                  //. timesrules
                  /. patternrules
                  /. w -> bunpairs];
  ];

(************************ CovDExpand **************************************)

defineTermwise[CovDExpand];

CovDExpand[ exp_ ] := TensorExpand[ exp //.
  Tensor[ t_, {i___}, {j___,U[k_],l___} ] :> 
    (Tensor[t,{i},{j,U[k],l}] /. LowerIndicesRule[U[k]]) //.
  Tensor[ t_, {i___}, {j___,k_} ] :> 
    expandedCovD[ Tensor[t,{i},{j}], {k} ]
  ];

expandedCovD[ x_, {k_,l__} ] := 
  expandedCovD[ expandedCovD[ x, {k} ], {l} ];

expandedCovD[ x_, {} ] := x;

expandedCovD[ x_, {k_} ] := Module[ {n,freeindices},
  freeindices = GetFreeIndices[x];
  Del[ Basis[k], x ] +
  Sum[ 
    (x /. freeindices[[n]] -> #) *
    If[IndexAltitude[#]===Upper,
      (*then contravariant index*)
          Tensor[Conn, { Pair[#], freeindices[[n]], k }, {} ],
      (*else covariant*)
        - Tensor[Conn, { freeindices[[n]], Pair[#], k }, {} ]
      ] & @ 
    dumdex[ freeindices[[n]] ],
    { n, Length[freeindices] }
    ]
  ];

(************************ ProductExpand ***********************************)

defineTermwise[ProductExpand];

ProductExpand[expr_] := 
  TensorExpand[ 
    expr //. 
      {HoldPattern[ Times[ x:( (_?(Rank[#]===1&)) |
                           ( (_?(Rank[#]===1&)) ^ (p_Integer) ) 
                         ).. ] ] :> 
         Module[{list},
           list = {x} //. {a___, b_ ^ (q_Integer?Positive), c___} :> 
                          Join[ {a}, Table[b,{q}], {c} ];
           (1/Length[list]!) *
           Plus @@ (TensorProduct @@ # &) /@ 
                   (list[[#]]&) /@ 
                   Permutations[Range[Length[list]]]
           ],
       HoldPattern[ Wedge[ x:(_?(Rank[#]===1&)).. ] ] :>
         (WedgeFactor[Length[{x}]] / Length[{x}]!) * Signature[{x}] *
         Plus @@ (Signature[#] * TensorProduct @@ # &) /@ 
                 ({x}[[#]] &) /@ 
                 Permutations[Range[Length[{x}]]]
      } //.
      {HoldPattern[ (x_ ? (Rank[#]===1&) ) ^ p_Integer?Positive ] :>
         TensorProduct @@ Table[ x, {p} ]}
  ];

(************************ CommuteCovD *************************************)

(* First check the arguments *)

CommuteCovD[] :=    Message[CommuteCovD::argrx, "CommuteCovD", 0, 3];
CommuteCovD[_] :=   Message[CommuteCovD::argr,  "CommuteCovD", 3];
CommuteCovD[_,_] := Message[CommuteCovD::argrx, "CommuteCovD", 2, 3];
CommuteCovD[_,_,_,extras__] := 
  Message[CommuteCovD::argrx, "CommuteCovD", 3+Length[extras], 3];

CommuteCovD[expr_,i_,j_] := 
   If[ !IndexQ[i,j],
       (*then*) Print["ERROR: invalid indices ",i,",",j,
                       " given to CommuteCovD"];
                Return[ERROR[CommuteCovD[expr,i,j]]],
       (*else*) Return[expr //. 
                  {Tensor[t_,{x___},{y___,i,j,z___}]
                   :> commutecd[Tensor[t,{},{}],{x},{y},{z},i,j]}]
     ];

(* commutecd does the actual computation.  *)

commutecd[Tensor[t_,{},{}], {a___}, {b___}, {c___}, i_, j_] :=
  Module[{n},
    (* First the curvature terms associated to undifferentiated indices *)
    Sum[ CovD[
           (Tensor[t,Join[Take[{a},n-1],{#},Drop[{a},n]],{b}]*
            InsertIndices[Curv, {{a}[[n]],Pair[#],i,j}]&) @
           dumdex[ {a}[[n]] ],
           {c} ],
         {n,Length[{a}]}] +
    (* Then curvature terms associated to previous derivative indices *)
    Sum[ CovD[
           (Tensor[t,{a},Join[Take[{b},n-1],{#},Drop[{b},n]]]*
            InsertIndices[Curv, {{b}[[n]],Pair[#],i,j}]&) @
           dumdex[ {b}[[n]] ], 
           {c} ],
         {n,Length[{b}]}] +
    (* Next the torsion terms *)
    Plus @@ (
      (CovD[#,{c}]&) /@ ( 
        (Tensor[t,{a},{b,#}] * InsertIndices[Tor, {Pair[#],i,j}]&) /@ 
        tandumdex[i]
        )
      ) +
    (* And finally the tensor itself with derivatives commuted. *)
    Tensor[t,{a},{b,j,i,c}]
  ];

(*********************** OrderCovD **************************************)

defineTermwise[OrderCovD];

OrderCovD[expr_] :=
  expr //.
  {Tensor[t_,{a___},{b___,i_,j_,c___}] :>
   commutecd[Tensor[t,{},{}],{a},{b},{c},i,j] /; !IndexOrderedQ[{i,j}]};

(********************* End of private context ***********************)

End[(*"`Private`"*)];

(********************************************************************)

(** Protect all the gadgets we've defined to guard against accidental
    modifications. **)

Protect["Ricci`*"];

(** We also unprotect all the global user parameters, which start with $ **)

Unprotect["Ricci`$*"];

EndPackage[ (* Ricci` *) ];
