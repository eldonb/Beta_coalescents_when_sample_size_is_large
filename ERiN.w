\pdfoutput=1
\documentclass[a4paper,10pt]{cweb}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
%%\usepackage[lf]{Baskervaldx}
%\usepackage[bigdelims,vvarbb]{newtxmath}
\usepackage{amsfonts, amsmath, amssymb}
\usepackage{fullpage}
\usepackage{marvosym}
\usepackage{bm}
\usepackage[round,numbers,super]{natbib}
\usepackage{color}
\usepackage{a4wide,fullpage}
\usepackage{setspace}
\usepackage{hyperref}
\hypersetup{
    colorlinks,
    linkcolor={gray!50!black},
    citecolor={blue!50!black},
    urlcolor={blue!80!black}
}
\usepackage{enumerate}
\usepackage{dsfont}
\usepackage[right]{lineno}
\usepackage{verbatim}
\usepackage{tabto}
\usepackage{lipsum}
\usepackage{orcidlink}
\setstretch{1.1}
\newcommand{\one}[1]{\ensuremath{\mathds{1}_{\left\{ #1 \right\}}}}
\newcommand{\EE}[1]{\ensuremath{\mathds{E}\left[ #1 \right]}}%
\newcommand{\im}{\ensuremath{\imath} }%
\newcommand{\jm}{\ensuremath{\jmath} }%
\newcommand{\prb}[1]{\ensuremath{\mathds{P}\left( #1 \right) } }%
\newcommand{\h}[1]{\ensuremath{\uptheta_{ #1 } } }%
\newcommand{\VV}[1]{\ensuremath{ \mathbb{V}\left( #1 \right)}}%
\newcommand{\hp}{\ensuremath{\theta_1}}
\newcommand{\hs}{\ensuremath{\theta_2}}
\newcommand{\D}{\ensuremath{\mathbb{D}}}
\newcommand{\F}{\ensuremath{\mathbb{F}} }
\newcommand{\G}{\ensuremath{\mathbb{G}} }
\newcommand{\bt}[1]{\textcolor{blue}{\tt #1}}
\newcommand{\aths}[1]{\textcolor{violet}{\textbf{\small #1}}}
\makeatletter
\renewcommand{\maketitle}{\bgroup\setlength{\parindent}{0pt}
\begin{flushleft}
  \textsf{\textbf{\@@title}}
\end{flushleft}
\begin{center}
  \textsc{\@@author}
\end{center}
\egroup
}
\makeatother
\title{Beta-coalescents when sample size is large \\ --- approximating $\EE{R_{i}^{N}(n)}$ }
\author{Bjarki Eldon\footnote{Email: \href{mailto:beldon11@@gmail.com}{beldon11@@gmail.com} \\ Supported by 
  Deutsche Forschungsgemeinschaft (DFG) - Projektnummer 273887127 
%% „funded by the Deutsche Forschungsgemein-schaft (DFG, German Research Foundation) –Projektnummer(n)“.
%% 
through DFG SPP 1819: Rapid Evolutionary Adaptation grant STE 325/17
to Wolfgang Stephan; acknowledge funding by the Icelandic Centre of
Research (Rann\'is)  through an Icelandic Research Fund (Ranns\'oknasj\'o{\dh}ur) Grant of Excellence no.\
185151-051 to Einar \'Arnason, Katr\'in Halld\'orsd\'ottir, Alison M.\
Etheridge, WS, and BE;  acknowledge  Start-up module grants
through SPP 1819 with Jere Koskela and Maite Wilke-Berenguer, and with
Iulia Dahmer. \\ \today}
\orcidlink{0000-0001-9354-2391} }

\begin{document}
\maketitle
\renewcommand{\abstractname}{\vspace{-\baselineskip}}

%%\rule{\textwidth}{.8pt}


\begin{abstract}
 Let  
$L_{i}^{N}(n) \equiv \sum_{j=1}^{\tau^{N}(n) } \# \left\{ \xi \in
\xi^{n,N}(j) : \#\xi = i \right\} $ and
$L^{N}(n) \equiv \sum_{j=1}^{\tau^{N}(n)} \# \xi^{N,n}(j) $ and
$\tau^{N}(n) \equiv \inf \left\{ j \in \mathds N : \# \xi^{n,N}(j) = 1
\right\} $ for $i \in \{1, 2, \ldots, n-1\}$, and
$R_{i}^{N}(n) \equiv L_{i}^{N}(n)/L^{N}(n) $ for $i=1,2,\ldots, n-1$
and $L^{N}(n) = \sum_{j}L_{j}^{N}(n) $. 
With this C++ code one estimates the functionals $\EE{R_{i}^{N}(n)}$
for $i = 1,2, \ldots, n-1$ of the ancestral process
$\left\{ \xi^{n,N}(r) : r \in \mathds N_{0} \right\}$ tracking the
random relations of $n$   sampled gene copies when the sample is from a
finite haploid panmictic population of constant size evolving
according to a given model of sweepstakes reproduction  (skewed
offspring number distribution).


\end{abstract}

\tableofcontents


@* {\bf Copyright}. 


Copyright {\textcopyright} {\the\year}  Bjarki Eldon \newline

This document and any source code it contains  is distributed under the terms of the GNU General Public Licence (version $\ge 3$).  You
should have received a copy of the licence along with this file (see file COPYING).  


    The source codes  described in this document  are  free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This document and the code it contains   is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this file (see COPYING).  If not, see \url{http://www.gnu.org/licenses/}.


@* {\bf Compilation,  output, and execution}. 
\label{compile}

 This CWEB
      \citep{knuth1994cweb} document (the {\tt .w} file) can be
      compiled with {\tt cweave} to generate a {\tt .tex} file, and
      with {\tt ctangle} to generate a {\tt .c} \citep{kernighan1988c}
      file.

One can use {\tt cweave} to generate a {\tt .tex} file, and {\tt
ctangle} to generate a {\tt .c} file. To compile the C++ code (the {\tt
.c} file), one needs the GNU Scientific Library.


Compiles on {\tt Linux Debian 6.12.9} with  {\tt ctangle 4.11} and {\tt g++ 14.2} and {\tt GSL 2.8}



Using a Makefile can be helpful, naming this file {\tt iguana.w}


 {\tt
iguana.pdf : iguana.tex \\
\tab\quad\quad\quad\quad cweave iguana.w \\
\tab\quad\quad\quad\quad        pdflatex iguana \\
\tab\quad\quad\quad\quad        bibtex iguana \\
\tab\quad\quad\quad\quad        pdflatex iguana \\
\tab\quad\quad\quad\quad        pdflatex iguana \\
\tab\quad\quad\quad\quad        ctangle iguana \\
\tab\quad\quad\quad\quad        g++ -Wall -Wextra -pedantic -std=c++26 -O3 -march=native -m64 -x c++ iguana.c -lm -lgsl -lgslcblas \\
        
       
clean :  \\
\tab\quad\quad\quad\quad        rm -vf iguana.c iguana.tex \\
}


Use {\tt valgrind} to check for memory leaks:

{\tt valgrind -v --leak-check=full --show-leak-kinds=all <program call>}


Use {\tt cppcheck} to check the code


{\tt cppchek  ---enable=all ----language=c++ <prefix>.c}

To generate estimates on a computer with several CPUs it may be
convenient to put  in a text file ({\tt simfile}):


{\tt ./a.out \$(shuf -i 484433-83230401 -n1) > resout<i>}


for  $i = 1,\ldots, y$ and use  {\tt
parallel}\cite{tange11:_gnu_paral}

{\tt parallel ---gnu -jy :::: ./simfile}



@* {\bf introduction}. 
\label{intro}


We consider a haploid population of fixed size $N$. Let
$X, X_1, \ldots, X_N$ be i.i.d.\ discrete random variables
taking values in $\{1, \ldots, \zeta(N) \}$; the $X_1, \ldots, X_N$
denote the random number of potential offspring  independently produced in a
given generation according to
\begin{equation}
\label{eq:1}
   \prb{X = k} = \frac{ (\zeta(N) +1)^\alpha }{ (\zeta(N)  + 1)^\alpha -
1 } \left( \frac{1}{k^\alpha} - \frac{1}{(k+1)^{\alpha}} \right),
\quad 1 \leq k \leq \zeta(N).
\end{equation}
The mass in \eqref{eq:1} is normalised so that
$\prb{ 1 \leq X \leq \zeta(N) } =1 $, and
$\prb{X = k} \ge \prb{X = k+1}$. Given a pool of at least $N$
potential offspring, we sample $N$ of them for the next generation
uniformly at random and without replacement.  Leaving out an atom at
zero gives $X_1 + \cdots + X_N \ge N$ almost surely, guaranteeing that
we always have at least $N$ potential offspring to choose from in each
generation.  If $1 < \alpha < 2$ and
$\liminf_{N\to \infty}\zeta(N)/N > 0$ the ancestral process tracing
the random ancestral relations of leaves converges in
finite-dimensional distributions to the
Beta$(\gamma, 2-\alpha,\alpha)$-coalescent with $0 < \gamma \le 1$; if
$\alpha \ge 2$ or $\zeta(N)/N \to 0$ (the ancestral process) converges
(in finite-dimensional distributions) to Kingman. Thus, the model
described in \eqref{eq:1} is a mathematically tractable model of
sweepstakes reproduction  (skewed offspring number distribution).



Let $\left\{\xi^{n,N}(j) : j\in \mathds N _{0} \right\} $ where
$\mathds N_{0} \equiv \{0,1,2, \ldots \}$ be the ancestral process
tracking the random ancestral relations of sampled gene copies when
the sample comes from a finite haploid panmictic population evolving
according to \eqref{eq:1}.  Let
$L_{i}^{N}(n) \equiv \sum_{j=1}^{\tau^{N}(n) } \# \left\{ \xi \in
\xi^{n,N}(j) : \#\xi = i \right\} $ and
$L^{N}(n) \equiv \sum_{j=1}^{\tau^{N}(n)} \# \xi^{N,n}(j) $ and
$\tau^{N}(n) \equiv \inf \left\{ j \in \mathds N : \# \xi^{n,N}(j) = 1
\right\} $ for $i \in \{1, 2, \ldots, n-1\}$, with $n$ being the
sample size.  Then
$L^{N}(n) = L_{1}^{N}(n) + \cdots + L_{n-1}^{N}(n)$; define
$R_{i}^{N}(n) \equiv L_{i}^{N}(n)/L^{N}(n)$.  We are interested in
$\EE{R_{i}^{N}(n)}$, and how $\EE{R_{i}^{N}(n)}$ behaves as a function
of $n$.  This has been investigated in the case of evolution according
to the Wright-Fisher model. We are interested in this question in the
case of evolution according to sweepstakes reproduction, with
numbers of potential offspring produced according to \eqref{eq:1} and
then sampled uniformly and without replacement.


The algorithm is summarised in \S~\ref{sec:code}, the code follows in
\S~\ref{sec:includes}--\S~\ref{sec:main}; we conclude in
\S~\ref{sec:concl}. Comments within the code are in \aths{this font and colour}

@* {\bf Code}. 
\label{sec:code}

The included libraries are listed in \S~\ref{sec:includes}, and  the
global constants in \S~\ref{sec:constants}.  The random number
generators are defined in \S~\ref{SEC:rng}.  The mass function in
Eq~\eqref{eq:1} is computed in \S~\ref{sec:mass},  the corresponding
CDF for sampling numbers of potential offspring is computed in \S~\ref{sec:cdf}.  A random
number of potential offspring  for a single individual is drawn in
\S~\ref{sec:randomjuvs}, and for all the $N$ individuals in
\S~\ref{sec:pool}. From the pool of potential offspring  (there are at least  $N$ of them
almost surely)   $N$ are sampled without replacement, and this is the
algorithm's 
Acciles's heel.  Sampling without replacement  is  inefficient.
To try to improve the efficiency in \S~\ref{sec:estimatecoalpr}  we  estimate the pairwise  coalescence
probability $c_{N}$ Eq~\eqref{eq:2}.  When there are two blocks
left, the random  time until the last two blocks merge is geometric with
success probability $c_{N}$.   When there are  three blocks left, we
can  use a similar approach,  i.e.\ estimate the corresponding
coalescence probabilities Eq~\eqref{eq:3} and Eq~\eqref{eq:4}, and
sample geometric times    In \S~\ref{sec:rmvhyper} we assign blocks
to  families given a realisation of number of potential offspring  per
individual.  If there are currently $n$ blocks, the  joint
distribution of number of blocks per family is multivariate
hypergeometric conditional on the number of potential offspring.  Given a
realisation $x_{1}\ldots, x_{N}$ of $X_{1}, \ldots, X_{N}$  we approximate
the  number of blocks per family with
\begin{equation}
\label{eq:5}
\prb{\nu = (v_{1}, \ldots, v_{N})} =  \frac{ \binom{x_{1}}{v_{1}} \cdots \binom{x_{N}}{v_{N}}    }{ \binom{x_{1} + \cdots + x_{N}}{n} } 
\end{equation}
 This approximation well approximates  first sampling $N$ potential offspring 
 and using the surviving  offspring in the hypergeometric
 Eq~\eqref{eq:5}.    The hypergeometric step in the algorithm is the
 bottleneck.  Given the number of blocks per family we can update the
 tree \S~\ref{sec:updatetree}. The tree is a vector of block sizes,
 and we merge blocks by first shuffling the order of the blocks, and
 then merging the rightmost blocks.  The new blocks are then
 appended to the  tree.  Thus, we only store the current configuration
 of the tree. Obviously if at most one block is assigned to each
 family then the tree is unchanged over the generation.    Given the
 current tree, the branch lengths for the current realisation  are updated
 in \S~\ref{sec:updateb}, and after each realisation the estimate of
 $\EE{R_{i}^{N}(n)}$ is updated in \S~\ref{sec:estimateri}.  We use
 the coalescence probability estimates \S~\ref{sec:estimatecoalpr} in
 \S~\ref{sec:threetwo}, in the case when there are at most three
 blocks left in the tree.  



Write $[n] = \{1,2, \ldots, n\}$ for any natural number $n$.  Let $n$
denote the sample size and $m \in [n]$ the current number of blocks.
We are  interested in the branch lengths and  require the
block sizes $(b_{1}, \ldots, b_{m})$ where $b_{j} \in [n]$ and
$b_{1} + \cdots + b_{m} = n$ are the current block  sizes

\begin{enumerate}
\item $\left( r_{1}(n), \ldots, r_{n-1}(n)\right) \leftarrow (0,\ldots, 0)$
\item for each of  $M$ experiments : 
\begin{enumerate}
\item set $(b_{1}, \ldots, b_{n}) \leftarrow  (1,\ldots, 1)$.  
\item set current branch lengths $\ell_{i}(n) \leftarrow  0$ for $i \in [n-1]$
\item set the current number of blocks  $m\leftarrow n$ 
\item {\bf while} $m > 1$ (at least two blocks; see  \S~\ref{sec:estimate}) 
\begin{enumerate}
\item update the current  branch lengths   $\ell_{b}(n) \leftarrow  1 +  \ell_{b}(n)$ for $b = b_{1},\ldots, b_{m}$;  \S~\ref{sec:updateb}
\item sample random number of  potential offspring $X_{1}, \ldots, X_{N}$ \S~\ref{sec:pool}
\item given $X_{1},\ldots, X_{N}$ assign blocks to families   \S~\ref{sec:rmvhyper}
\item merge blocks assigned to the same family  and update current number of blocks \S~\ref{sec:updatetree}
\end{enumerate}
\item given a realisation $\ell_{1}(n), \ldots, \ell_{n-1}(n)$   of branch lengths update the estimate of
$\EE{R_{i}^{N}(n)}$  \S~\ref{sec:estimateri} ; $r_{i}(n) \leftarrow r_{i}(n) + \ell_{i}(n)/\sum_{j}\ell_{j}(n)$   
\end{enumerate}
\item return an estimate  $(1/M)r_{i}(n)$  of   $\EE{R_{i}^{N}(n)}$ for $i = 1,2, \ldots, n-1$
\end{enumerate}




@*1 {\bf includes}.
\label{sec:includes}

the included libraries; we use the GSL Library


@<includes@>=@#
#include <iostream>
#include <vector>
#include <random>
#include <functional>
#include <memory>
#include <utility>
#include <algorithm>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <list>
#include <string>
#include <fstream>
#include <forward_list>
#include <chrono>
#include <assert.h>
#include <math.h>
#include <unistd.h>
#include <omp.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_gamma.h>


@*1 {\bf constants}.
\label{sec:constants}


the global constants

@<constants@>=@#

/* \newline  \aths{the $\alpha$ parameter in \eqref{eq:1} }  */
const double CONST_ALPHA = 1.0 ;
/* \newline \aths{ population size $N$}  */
const size_t CONST_POP_SIZE = 1e3 ;
/* \newline \aths{ the cutoff $\zeta(N)$ in \eqref{eq:1}}  */
const double CONST_CUTOFF = 1.0e3  ;
/* \newline \aths{ sample size} */
const size_t CONST_SAMPLE_SIZE = 10 ;
/* \newline \aths{ number of experiments}  */
const double CONST_NUMBER_EXPERIMENTS = 2500. ;


@*1 {\bf the random number generator}. 
\label{SEC:rng}


@<gslrng@>=@#
/* \newline \aths{ the GSL random number engine }  */
gsl_rng * rngtype ;
 /* \newline \aths{ obtain a seed out of thin air for the random number engine}  */
 std::random_device randomseed;
  /* \newline \aths{ Standard mersenne twister  random number engine seeded with |randomseed()|} */
  std::mt19937_64 rng(randomseed());

/* \newline \aths{ set up and initialise the GSL random number generator}  */
static void setup_rng( unsigned long int s )
{
        const gsl_rng_type *T ; 
        gsl_rng_env_setup(); 
        T = gsl_rng_default ;
        rngtype = gsl_rng_alloc(T);
        gsl_rng_set( rngtype,  s) ;
}



@*1 {\bf the mass function}. 
\label{sec:mass}


The mass function in \eqref{eq:1}

@<mass@>=@#
static double massfunction( const double j)
{

@#

 /* \newline  \aths{ |CONST_CUTOFF|  and |CONST_ALPHA|   \S~\ref{sec:constants}; recall $\zeta(N)$ from   \eqref{eq:1} } */
return ( (pow(1.0/j, CONST_ALPHA) -   pow(1.0/(1. + j),
CONST_ALPHA))/( 1. -  pow( 1./(CONST_CUTOFF + 1.), CONST_ALPHA)) ) ;
}


@*1 {\bf cdf}.
\label{sec:cdf}

 define the function for computing the 
   CDF  for the  distribution of number of potential offspring  in \eqref{eq:1}.

@<cdf@>=@#
static void mass_function( std::vector<double>& vcdf )
{

@#

  vcdf.clear() ;
  /* \newline \aths{ we take $\prb{X = 0} = 0$} */
  vcdf.push_back(0.0) ;
  /* \newline  \aths{ |CONST_CUTOFF|   \S~\ref{sec:constants}; recall $\zeta(N)$ from   \eqref{eq:1} } */
  for( double j = 1; j <= CONST_CUTOFF ; ++j){
  /* \newline \aths{ adding upp the mass function \S~\ref{sec:mass}} */
    vcdf.push_back( vcdf.back() + massfunction(j) );}
}


@*1 {\bf random number of potential offspring}.
\label{sec:randomjuvs}

the function for sampling a random number of potential offspring
returning $\min \left\{ j \in \mathds{N} : F(j) \ge u \right  \} $ for $u$ a random
uniform, and $F$ the CDF computed in \S~\ref{sec:cdf}

@<randomjuvs@>=@#
static size_t sample_juveniles( const std::vector<double>& vcdf  )
{
@#

  size_t  j = 1;
  const double u = gsl_rng_uniform( rngtype ); 
  while( vcdf[j] < u){
    ++j ; }
  return(j) ;
}

@*1 {\bf random unbounded potential offspring}.
\label{sec:unboundedx}

unbounded; this is here for completeness more than anything else, one
could interpret $\zeta(N) = N\log N$ as ``unbounded'' since then
$\zeta(N)/N\to \infty$ as $N\to \infty$

@<unbounded@>=@#
static std::size_t randomX()
{

/* \newline \aths{ |CONST_ALPHA| from \S~\ref{sec:constants}} */

   return floor( 1./pow( gsl_rng_uniform_pos(rngtype),
   1./CONST_ALPHA)) ;
}

@*1 {\bf pool of potential offspring}. 
\label{sec:pool}

sample a random pool of potential offspring, i.e.\ each of the $N$ individuals
independently contributes a random number of potential offspring  according to
\eqref{eq:1} 

@<pool@>=@#
static size_t sample_pool_juveniles( std::vector<size_t>& pool_juvs, const std::vector<double>& v_cdf)
{

@#

  pool_juvs.clear();
  size_t s = 0;
  for( size_t i = 0; i < CONST_POP_SIZE; ++i){
  /* \newline  \aths{ record the random number of potential offspring  for each individual \S~\ref{sec:randomjuvs}}   */
  /* \newline \aths{  |sample_juveniles(v_cdf)| for bounded} */
  /* \newline  \aths{ |randomX()| for  unbounded} */
    pool_juvs.push_back( sample_juveniles(v_cdf) );
    s += pool_juvs.back() ; }
  /* \newline  \aths{  return  the total number of juveniles $X_{1} + \cdots + X_{N}$ } */
  return (s);
}



@*1 {\bf estimate coalescence probabilities}. 
\label{sec:estimatecoalpr}


estimate coalescence probabilities for  speeding up reaching the most
recent common ancestor.  When only two blocks left we can sample a
geometric with success probability the pairwise coalescence
probability.  When only three blocks left can sample between a
pairwise merger and  a triple merger. Given a realisation $x_{1},
\ldots, x_{N}$ of $X_{1},
\ldots, X_{N}$ with $s_{N} := x_{1} + \cdots + x_{N}$  the pairwise coalescence
probability is
\begin{equation}
\label{eq:2}
c_{N} =  \sum_{j=1}^{N} \frac{x_{j}(x_{j} - 1)}{s_{N}(s_{N}-1)},
\end{equation}
a 3-merger when three blocks is
\begin{equation}
\label{eq:3}
c_{N}(3;3) = \sum_{j=1}^{N} \frac{(x_{j})_{3}}{(s_{N})_{3}},
\end{equation}
a 2-merger when three blocks is
\begin{equation}
\label{eq:4}
c_{N}(3;2) =  \sum_{j=1}^{N} \frac{3(x_{j})_{2}(s_{N} - x_{j})}{(s_{N})_{3}}.
\end{equation}

@<estimatecoalpr@>=@#
static void estimate_coalescence_probabilities( std::vector<double>& v_cN, const std::vector<double>& v_cdf,   std::vector<size_t>& v_pool_jvs)
{

@#

  size_t SN {} ;
  /* \newline \aths{ estimate the coalescence probabilites from $10^{3}$ experiments}  */
  for( size_t i = 0 ; i < 1000 ; ++i){
  /* \newline \aths{ sample a pool of potential offspring and record the total number of
  them  \S~\ref{sec:pool}}  */
    SN = sample_pool_juveniles( v_pool_jvs,  v_cdf) ;
    for( size_t j = 0 ; j < CONST_POP_SIZE; ++j){
      /* \newline \aths{ the pairwise probability \eqref{eq:2}} */
      v_cN[0] += (static_cast<double>( v_pool_jvs[j])/static_cast<double>(SN)) * (static_cast<double>(v_pool_jvs[j] - 1) / static_cast<double>(SN - 1)) ;
      /* \newline \aths{ a 3-merger when three blocks \eqref{eq:3}}  */
      v_cN[1] += (static_cast<double>( v_pool_jvs[j])/static_cast<double>(SN)) *  (static_cast<double>( v_pool_jvs[j]-1)/static_cast<double>(SN-1)) *  (static_cast<double>( v_pool_jvs[j]-2)/static_cast<double>(SN-2)) ;
      /* \newline \aths{ a merger of two of three blocks Eq~\eqref{eq:4}}   */
      v_cN[2] += 3.*(static_cast<double>(SN - v_pool_jvs[j]) / static_cast<double>(SN)) * (static_cast<double>( v_pool_jvs[j])/static_cast<double>(SN - 1)) *(static_cast<double>(v_pool_jvs[j] - 1) / static_cast<double>(SN-2)) ;
    } }
  /* \newline \aths{ take the average } */
  v_cN[0] /= 1000. ;
  v_cN[1] /= 1000. ;
  v_cN[2] /= 1000. ;

}


@*1 {\bf sample a multivariate hypergeometric }.
\label{sec:rmvhyper}

assign the ancestral blocks to families; given a realisation of random
numbers of potential offspring  the joint  number of blocks per family is a
multivariate hypergeometric.  We sample the marginals and update. 

@<rmvhyper@>=@#
static void rmvhyper( std::vector<size_t>& merger_sizes,  size_t k, const std::vector<size_t>& v_juvs, const size_t SN, gsl_rng *r )
{

@#


  /* \newline \aths{ |k| is the current number of lines}  */
  merger_sizes.clear();
  size_t number_of_new_lines = 0 ;
  size_t n_others =  SN - v_juvs[0] ;
  /* \newline \aths{ sample the number of blocks assigned to the first family }  */
  size_t x = gsl_ran_hypergeometric( r, v_juvs[0], n_others, k);
  if( x > 1){
    /* \newline \aths{ only record  merger sizes} */
    merger_sizes.push_back(x ); }
    /* \newline \aths{ update the remaining number of blocks}  */
  k -= x ;
  /* \newline \aths{ update new number of lines}  */
  number_of_new_lines += ( x > 0 ? 1 : 0) ;
  size_t i =0 ;
    /* \newline \aths{ we can stop as soon as all lines 
       have been assigned to a family } */
  while( (k > 0) && (i < CONST_POP_SIZE-1) ){
  /* \newline \aths{ set the index to the one being sampled from}  */
    ++i ;
    /* \newline \aths{ update |n_others| }  */
    n_others -= v_juvs[i] ;
    x =  gsl_ran_hypergeometric( r,  v_juvs[i], n_others, k );
    if( x > 1){
      merger_sizes.push_back( x) ; }
      /* \newline  \aths{ update the remaining number of blocks} */
      k -= x ;
    /* \newline \aths{ update new number of lines}  */
    number_of_new_lines += ( x > 0 ? 1 : 0) ;
  }
  /* \newline \aths{ check if at least two lines assigned to last individual} */
  if( k > 1){
    merger_sizes.push_back( k); }
  /* \newline \aths{ check if at least one line assigned to last individual}  */
  if( k > 0){
    number_of_new_lines += 1;}
}


@*1 {\bf update the tree}. 
\label{sec:updatetree}

update the tree; the tree is a vector of block sizes. If there are
mergers we shuffle the tree and then   consequtively  merge blocks by summing and recoding  the size of the
merging blocks in each merger, removing the blocks that merge and
eventually adding the new blocks to the tree, the rightmost blocks
merging each time.  

@<updatetree@>=@#
static void update_tree( std::vector<size_t>& tree, const std::vector<size_t>& merger_sizes )
{

@#

  std::vector<size_t> new_blocks {} ;
  if( merger_sizes.size() > 0){
    /* \newline \aths{ at least one merger } */
    new_blocks.clear() ;
    /* \newline \aths{ shuffle the tree}  */
    std::ranges::shuffle( tree,  rng );
    /* \newline \aths{ loop over the mergers} */
    for( const auto &m: merger_sizes){
      /* \newline \aths{  |m| is number of blocks merging;  |m| is at least two; append
      new block to vector of  new blocks}  */
      assert( m > 1) ;
      /* \newline \aths{ record the size of the new block by summing the sizes of the
      merging blocks}  */
      new_blocks.push_back( std::accumulate( std::rbegin( tree), std::rbegin(tree) + m, 0) ) ;
      assert( new_blocks.back() > 1) ;
      /* \newline \aths{ remove the rightmost |m| merged  blocks from tree } */
      tree.resize( tree.size() - m) ;
    }
    /* \newline \aths{ append new blocks to tree} */
    tree.insert( tree.end(),  new_blocks.begin(),  new_blocks.end() ) ;
  }
  /* \newline \aths{ if no mergers then tree is unchanged} */
}


@*1 {\bf update the branch lengths}. 
\label{sec:updateb}

update the branch lengths

@<updateb@>=@#
static void update_ebib( const std::vector< size_t>& tree,  std::vector<double>& vebib)
{

@#

  for( const auto &b: tree){
    /* \newline \aths{ |b| is size of current block} */
    /* \newline \aths{ update the total tree size and then the branch length
    corresponding to the size of the block} */
    vebib[0] += 1.0 ;
    vebib[ b ]  += 1.0 ;}
}




@*1 {\bf update estimate of $\EE{R_{i}^{N}(n)}$}.
\label{sec:estimateri}

update estimate of  $\EE{R_{i}^{N}(n)}$ given branch lengths from one
realisation of a tree 

@<updateri@>=@#
static void update_estimate_ebib( const std::vector<double>& v_tmp,  std::vector<double>& v_ebib)
{

@#


  for ( size_t i = 1 ; i < CONST_SAMPLE_SIZE ; ++i){
    v_ebib[i] += v_tmp[i]/v_tmp[0] ;}
}



@*1 {\bf three or two blocks left}. 
\label{sec:threetwo}

at most three blocks left, so sample times  using the estimates of the
coalescence probabilities \S~\ref{sec:estimatecoalpr}


@<threetwo@>=@#
static void three_or_two_blocks_left(  std::vector<double>& tmp_bib, const std::vector<double>& v_cN,  std::vector<size_t>& v_tree)
{

  double Tk = 0.;
  double Tkk = 0. ;
  size_t newblock {} ;
  switch( v_tree.size() ){
  case 3 : {
    /* \newline \aths{ three lines left so  sample the two waiting times for a 3-merger and a 2-merger} */
    Tk = static_cast<double>( gsl_ran_geometric(rngtype, v_cN[1] ) ) ;
    Tkk = static_cast<double>( gsl_ran_geometric(rngtype, v_cN[2]));
    if( Tk < Tkk){
      /* \newline \aths{ all three blocks merge;  update the branch lengths} */

      tmp_bib[0] += (3. * Tk) ;
      tmp_bib[ v_tree[0]] += Tk ;
      tmp_bib[ v_tree[1]] += Tk ;
      tmp_bib[ v_tree[2]] += Tk ;
      /* clear the tree */
      v_tree.clear() ;
      assert( v_tree.size() < 1);
    }
    else{
      /* \newline \aths{ a 2-merger occurs followed by a merger of the last two
      blocks } */
      tmp_bib[0] += (3. * Tkk) ;
      tmp_bib[ v_tree[0]] += Tkk ;
      tmp_bib[ v_tree[1]] += Tkk ;
      tmp_bib[ v_tree[2]] += Tkk ;
      /* \newline \aths{ shuffle the tree } */
      std::ranges::shuffle( v_tree,  rng );
      newblock = v_tree[1] + v_tree[2] ;
      v_tree.resize(1) ;
      v_tree.push_back( newblock);
      assert( v_tree.size() == 2 );
      /* \newline \aths{ sample waiting time until merger of last two blocks} */
      Tk = static_cast<double>( gsl_ran_geometric(rngtype, v_cN[0] ) ) ;
      tmp_bib[0] += (2. * Tk) ;
      tmp_bib[ v_tree[0]] += Tk ;
      tmp_bib[ v_tree[1]] += Tk ;
      v_tree.clear() ;
      assert( v_tree.size() < 1);
    }
    break ; }
  case 2 : {
    /* \newline \aths{ two blocks left} */
  Tk = static_cast<double>( gsl_ran_geometric(rngtype, v_cN[0] ) ) ;
    tmp_bib[0] += (2. * Tk) ;
    tmp_bib[ v_tree[0]] += Tk ;
    tmp_bib[ v_tree[1]] += Tk ;
    v_tree.clear() ;
    assert( v_tree.size() < 1);
    break ; }
  default : break ;
  }
}

@*1 {\bf estimate $\EE{R_{i}^{N}(n)}$ }. 
\label{sec:estimate}

estimate   $\EE{R_{i}^{N}(n)}$ from a given number of experiments. 

@<theestimator@>=@#
static void estimate_ebib( )
{

@#

  std::vector< double > v_cdf ;

@#

  v_cdf.reserve( static_cast<size_t>(CONST_CUTOFF) + 1) ;
  /* \newline \aths{ compute the CDF function for sampling potential offspring  \S~\ref{sec:cdf}} */
  mass_function( v_cdf);

@#

  std::vector<size_t> v_number_juvs ;
  v_number_juvs.reserve( CONST_POP_SIZE);

@#

/* \newline \aths{ the tree  (current block sizes); initially all blocks are singletons (of size 1)} */
  std::vector<size_t> v_tree (CONST_SAMPLE_SIZE, 1) ;

@#

  std::vector<size_t> v_merger_sizes {};


@#

  v_merger_sizes.reserve( CONST_SAMPLE_SIZE );


@#

  std::vector<double> v_tmp_ebib (CONST_SAMPLE_SIZE, 0.0) ;

@#

  std::vector<double> v_ebib (CONST_SAMPLE_SIZE, 0.0) ;


@#

  std::vector<double> v_coal_probs (3, 0.0) ;

  /* \newline \aths{ estimate the coalescence probs \S~\ref{sec:estimatecoalpr}} */
  estimate_coalescence_probabilities( v_coal_probs, v_cdf, v_number_juvs) ;
  
  size_t SN = 0;
  double number_experiments = CONST_NUMBER_EXPERIMENTS + 1.;
  while( --number_experiments  > 0.){
  /* \newline \aths{ initialise the tree as all singletons} */
  v_tree.clear();
  v_tree.assign( CONST_SAMPLE_SIZE, 1);
  /* \newline \aths{ initialise the container for the branch length for  the current
  realisation} */
  std::fill(  std::begin( v_tmp_ebib), std::end(  v_tmp_ebib ), 0.0 );
  assert( std::accumulate(  std::begin( v_tmp_ebib), std::end(
  v_tmp_ebib ), 0.0) == 0); 
  while( v_tree.size() > 1){
    
    /* \newline \aths{ record the branch lengths for the current tree configuration} */
    update_ebib( v_tree, v_tmp_ebib) ;
    if( v_tree.size() > 3 ){
    /* \newline  \aths{ sample pool of potential offspring} */
    SN = sample_pool_juveniles( v_number_juvs, v_cdf) ;
    /* \newline \aths{ compute the merger sizes \S~\ref{sec:rmvhyper}} */
    rmvhyper( v_merger_sizes, v_tree.size(),  v_number_juvs, SN, rngtype) ;
    /* \newline \aths{ update the tree \S~\ref{sec:updatetree}}  */
    update_tree( v_tree, v_merger_sizes);}
    else{
      /* \newline \aths{ at most three blocks left \S~\ref{sec:threetwo}} */
      three_or_two_blocks_left( v_tmp_ebib, v_coal_probs, v_tree) ;
    }
  }
  /* \newline \aths{ update estimate of   $\EE{R_{i}^{N}(n)}$}  */
  update_estimate_ebib( v_tmp_ebib, v_ebib);
  }
  /* \newline \aths{ print the  estimate of $\EE{R_{i}^{N}}$} */
  for( const auto&r: v_ebib){
    std::cout << r << '\n' ; }
}
 


@*1 {\bf the main module}.
\label{sec:main}

The |main| function

@C

@<includes@>@#
@<gslrng@>@#
@<constants@>@#
@<mass@>@#
@<cdf@>@#
@<randomjuvs@>@#
@<unbounded@>@#
@<pool@>@#
@<estimatecoalpr@>@#
@<rmvhyper@>@#
@<updatetree@>@#
@<updateb@>@#
@<updateri@>@#
@<threetwo@>@#
@<theestimator@>@#

int main(int argc, char *argv[])
{

/* \newline \aths{ initialise the GSL random number generator |rngtype|
\S~\ref{SEC:rng} } */
setup_rng(  static_cast<unsigned long int>(atoi(argv[1])) );

/* \newline \aths{ estimate $\EE{R_{i}^{N}(n)}$ \S~\ref{sec:estimate}} */
estimate_ebib( ) ;

gsl_rng_free( rngtype ); 
return GSL_SUCCESS ; 
}


@* {\bf conclusion and bibliography}.
\label{sec:concl}


In the Schweinsberg model\cite{schweinsberg03} it is assumed that
individuals can produce arbitrarily many potential offspring. From the
pool of all potential offspring produced at any given time the
surviving offspring are sampled without replacement (provided there is
enough of them) to survive and replace the parents.  In \eqref{eq:1}
an upper bound $\zeta(N)$ is introduced to the distribution of number
of potential offspring. It can be shown that if $\zeta(N)/N \to K$
where $K >0 $ is fixed then the ancestral process converges to a
Beta$(\gamma,2-\alpha,\alpha)$-coalescent \cite{chetwyn-diggle_beta},
a truncated (or incomplete) form of the
Beta$(2-\alpha,\alpha)$-coalescent of \cite{schweinsberg03}.  With
this C++ code we estimate the functionals $\EE{R_{i}^{N}(n)}$;
estimates of $\EE{R_{i}^{N}(n)}$ may then be compared to estimates of
$\EE{R_{i}(n)}$, the functionals predicted by a given coalescent.





%%\bibliographystyle{plain}
%%\bibliography{refs}

\begin{thebibliography}{1}

\bibitem{chetwyn-diggle_beta}
JA~Chetwyn-Diggle, Bjarki Eldon, and Matthias Hammer.
\newblock Beta-coalescents when sample size is large.
\newblock In preparation, 2025+.

\bibitem{kernighan1988c}
Brian~W Kernighan and Dennis~M Ritchie.
\newblock The {C} programming language, 1988.

\bibitem{knuth1994cweb}
Donald~Ervin Knuth and Silvio Levy.
\newblock {\em The CWEB system of structured documentation: version 3.0}.
\newblock Addison-Wesley Longman Publishing Co., Inc., Reading, Massachusetts,
  1994.

\bibitem{schweinsberg03}
J~Schweinsberg.
\newblock Coalescent processes obtained from supercritical {G}alton-{W}atson
  processes.
\newblock {\em Stoch Proc Appl}, 106:107--139, 2003.

\bibitem{tange11:_gnu_paral}
O~Tange.
\newblock {GNU} parallel -- the command-line power tool.
\newblock The USENIX Magazine, 2011.

\end{thebibliography}






@
\end{document}
