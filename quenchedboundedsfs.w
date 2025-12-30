\pdfoutput=1
\documentclass[a4paper,10pt]{cweb}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
%%\usepackage[lf]{Baskervaldx}
%\usepackage[bigdelims,vvarbb]{newtxmath}
\usepackage{amsfonts, amsmath}
\usepackage{fullpage}
\usepackage{marvosym}
\usepackage{bm}
%%\usepackage[round,numbers,super]{natbib}
\usepackage{natbib}
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
\usepackage{url}
\usepackage{orcidlink}
\setstretch{1}
\usepackage[position=top,font={sf,md,up}]{subfig}
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
\newcommand{\NN}{\ensuremath{\mathds{N}} }
\newcommand{\set}[1]{\ensuremath{\left\{ #1 \right\}}}
\newcommand{\svigi}[1]{\ensuremath{\left( #1 \right)}}
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
} \makeatother
\title{Beta-coalescents when sample size is large \\ --- approximating $\EE{\widetilde R_{i}^{N}(n)}$}
\author{Bjarki Eldon\footnote{
\href{mailto:beldon11@@gmail.com}{beldon11@@gmail.com}}\footnote{Funding by Icelandic Centre of Research
(Rann\'is) through an Icelandic Research Fund
(Ranns\'oknasj\'o{\dh}ur) Grant of Excellence no.\ 185151-051 to Einar
\'Arnason, Katr\'in Halld\'orsd\'ottir, Wolfgang Stephan, Alison
Etheridge, and BE; DFG SPP 1819 Programme Rapid Evolutionary
Adaptation Start-up module grants with Jere Koskela, Maite Wilke
Berenguer,  and with Iulia Dahmer \\ \today }\orcidlink{0000-0001-9354-2391} }



\begin{document}
\maketitle
\renewcommand{\abstractname}{\vspace{-\baselineskip}}


%%\rule{\textwidth}{.8pt}


\begin{abstract}
With this C++ code one estimates mean relative branch lengths
$\EE{\widetilde R_{i}^{N}(n)}$ where $R_{i}^{N}(n) =
L_{i}^{N}(n)/\sum_{j=1}^{n-1}L_{j}^{N}(n)$ and $L_{i}^{N}(n)$ is the
random branch length supporting $i \in \set{1,2,\ldots, n-1}$ leaves
when the sample comes from a finite haploid panmictic population of
constant size evolving according to sweepstakes reproduction (skewed
offspring number distribution) and time is measured in discrete time
steps (generations).  The random sample of $n$ leaves is from a finite
haploid panmictic population of constant size $N$.  The population
evolves according to a model of sweepstakes reproduction.  The key
point is the conditioning on an ancestry (represented by the
sigma-field $\mathcal{A}^{(N,n)}$ recording the relations between
individuals); the population evolves forward in time and at each time
step (generation) a random sample is drawn from the population and the
sample tree checked for completeness; when a complete sample tree is
obtained, i.e.\ when a common ancestor for all the leaves in the
sample is found, the sample has a fixed ancestry, a fixed tree.  The
sample tree is traced and the branch lengths recorded. This process is
then repeated a given number of times, each time starting from scratch
with a new population, thus averaging over random ancestries (random
complete sample trees).  In this way one obtains an estimate of
$\EE{\widetilde R_{i}^{N}(n)}$ for $i = 1,2,\ldots, n-1$.
\end{abstract}


\tableofcontents


@* {\bf Copyright}.
\label{sec:copyright}

Copyright {\textcopyright} {\the\year}  Bjarki Eldon \newline

This document and any source code it contains is distributed under the
terms of the GNU General Public Licence (version $\ge 3$).  You should
have received a copy of the licence along with this file (see file
COPYING).


    The source codes described in this document are free software: you
    can redistribute it and/or modify it under the terms of the GNU
    General Public License as published by the Free Software
    Foundation, either version 3 of the License, or (at your option)
    any later version.

    This document and the code it contains is distributed in the hope
    that it will be useful, but WITHOUT ANY WARRANTY; without even the
    implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
    PURPOSE.  See the GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this file (see COPYING).  If not,
    see \url{http://www.gnu.org/licenses/}.


@* {\bf compilation and output}. 
\label{sec:compile}


 This CWEB
      \citep{knuth1994cweb} document (the {\tt .w} file) can be
      compiled with {\tt cweave} to generate a {\tt .tex} file, and
      with {\tt ctangle} to generate a {\tt .c} \citep{kernighan1988c}
      file.

One can use {\tt cweave} to generate a {\tt .tex} file, and {\tt
ctangle} to generate a {\tt .c} file. To compile the C++ code (the {\tt
.c} file), one needs the GNU Scientific Library.



Compiles on {\tt Linux Debian trixie/sid} with kernel {\tt 6.12.10-amd64}  and   {\tt ctangle 4.11} and {\tt g++ 14.2} and {\tt GSL 2.8}



{\tt  g++ -Wall -Wextra -pedantic -std=c++26 -O3 -march=native -m64 -x c++ <prefix>.c -lm -lgsl -lgslcblas }


Use {\tt valgrind} to check for memory leaks:

{\tt valgrind -v ---leak-check=full ---leak-resolution=high ---num-callers=40 ---vgdb=full <program call>}


Use {\tt cppcheck} to check the code:

{\tt cppchek  ---enable=all ----language=c++ <prefix>.c}


To generate estimates on a computer with several CPUs it may be
convenient to put  in a text file ({\tt simfile}):


{\tt ./a.out \$(shuf -i 484433-83230401 -n1) > resout<i>}


for  $i = 1,\ldots, y$ and use  {\tt
parallel}\cite{tange11:_gnu_paral}

{\tt parallel ---gnu -jy :::: ./simfile}




@* {\bf intro}.
\label{sec:intro}



A coalescent is a probabilistic description of the random ancestral
relations of sampled gene copies (leaves).  A coalescent
$\{\xi\} \equiv \{ \xi(t); t \ge 0\}$ is a Markov chain on the
partitions of $\mathds{N} = \{1,2,\ldots\} $, where the only
transitions are the merging of blocks (elements of a partition);
restricting to $n\in N$ gives a coalescent \set{\xi^{n}(t) : t \ge 0 }
on the partitions of $[n] = \{1,2,\ldots, n\}$.  We will consider
coalescents describing the ancestral relations of gene copies of a
single non-recombining locus (a contiguous non-recombining segment of
a chromosome) in a single haploid panmictic population of constant
size $N$.  A ``quenched coalescent'' would be a coalescent obtained by
conditioning on a random ancestry  of the individuals in
the population.


For the population model we consider an extension of the Schweinsberg model
\citep{schweinsberg03}. Let $X$ denote the random number of  potential
offspring of an arbitrary individual. Then, with $1 < \alpha < 2$,    
\begin{equation}
\label{eq:1}
\prb{X\ge k} = \frac{1}{k^{\alpha}}, \quad k \in \set{1,2,\ldots}
\end{equation}
is the unbounded distribution, and
\begin{equation}
\label{eq:2}
\prb{X=k} =  C\left( \frac{1}{k^{\alpha}} -  \frac{1}{(1+k)^{\alpha}}  \right), \quad k \in \set{1,2,\ldots, \zeta(N)}
\end{equation}
with $\zeta(N)$ the upper bound on the distribution, with $C$ so that
$\prb{1 \le X\le \zeta(N)} = 1$.  The population evolves according to sweepstakes reproduction. 
 In  each generation the current individuals independently
produce potential offspring according to \eqref{eq:1} (or
\eqref{eq:2}), from the pool of all potential offspring $N$ of them
are sampled uniformly and without replacement to replace the current
individuals (since $\prb{X \ge 1} = 1$ we will have $N$ potential
offspring each and every time almost surely).

Let $L_{i}^{N}(n)$ denote the random total  length of branches
supporting $i \in \set{1,2,\ldots, n-1}$ leaves; then $L^{N}(n) =
L_{1}^{N}(n) + \cdots + L_{n-1}^{N}(n)$  is the total tree size, and
$R_{i}^{N}(n) = L_{i}^{N}(n)/L^{N}(n)$ the relative branch length.
The quantity $R_{i}^{N}(n)$ is well defined since $L^{N}(n) \ge n$
almost surely.    We estimate $\EE{\widetilde R_{i}^{N}(n)}$
where  the sigma-field  $\mathcal{A}^{(N,n)}$ is the random ancestry
of all  the $N$  individuals in the population and information on which
$n$ individuals (leaves)  were sampled.  The coalescence probability
can be shown to be different  between the  annealed and the quenched
coalescent  when  the two gene copies  come from a population evolving
according to a specific model of  sweepstakes reproduction 
\citep{Diamantidis2024}.


In Figure~\ref{fig:quenchedannealedsfs} in \S~\ref{sec:examples} we
compare estimates of $\EE{\widetilde R_{i}^{N}(n)}$ to estimates of
$\EE{R_{i}(n) }$. It may not be a perfectly appropriate comparison but
at this time we do not have a quenched coalescent. The comparison can
still inform about how much $\EE{\widetilde R_{i}^{N}(n)}$ deviates
from $\EE{R_{i}(n) }$ as predicted by the annealed coalescent, e.g.\
as derived in \citep{schweinsberg03}.



The algorithm is summarised in \S~\ref{sec:code}, the code follows in
\S~\ref{sec:includes}--\S~\ref{sec:main}; we conclude in
\S~\ref{sec:concl}. Comments within the code are in \aths{this font and colour}


@* {\bf code}.
\label{sec:code}


The population tree is recorded as individuals living on levels, and
each individual ``points to'' the level of its immediate ancestor.  
\begin{center}
generation :  levels \\
0:   1    2     3     4     5     \textcolor{red}{6}    7 8    9    10 \newline
1:   \textcolor{red}{6}    6    9    8   \textcolor{red}{6}    6   10    1    1     \textcolor{blue}{2} \newline 
2:   4    \textcolor{red}{5}  10    9    6    \textcolor{red}{1}    5   \textcolor{blue}{10}   \textcolor{blue}{10}     2 \newline
\end{center}
Write $A_{\ell}(g)$ for the ancestor of individual on level $\ell$ in
generation $g$. One may imagine the individuals are assigned levels,
one individual on each level, and there are $N$ levels.  Suppose we
sample the individual on level 2 in generation 2, then $A_{2}(2) = 5$,
and $A_{5}(1) = 6$.  Suppose we have also sampled the individual on
level 6 in generation 2; then $A_{6}(2) = 1$, and $A_{1}(1) = 6$. The
two sampled individuals share the ancestor on level 6 in generation
0. The ancestry of the two individuals is shown in red in the tree
above for $N=10$. The ancestry of the individuals on levels 8 and 9 in
generation 2 is shown in blue. Each individual points  to
its immedate ancestor, and this is sufficient information to trace the
ancestry of any given individual and find a  common ancestor for any set of individuals.  


@*1 {\bf a summary of the algorithm}.
\label{sec:pseudocode}

In this section we summarize the algorithm for estimating
$\EE{\widetilde R_{i}^{N}(n)}$.  We repeatedly simulate a population
ancestry from scratch and regularly search for an ancestor of a random
sample of leaves.


The algorithm may well be improved; the focus here is not on
efficiency but correctness, to put together a working algorithm that
does the correct thing.

\begin{enumerate}
\item initialise the population ancestry  with the indexes $0,1,\ldots,
N-1$
\item draw uniformly without replacement  $n$ levels  from
$\set{0,1,2, \ldots,N-1}$  representing the
levels  of the leaves  of a new sample \S~\ref{sec:randomsample}
\item until a complete sample tree is found  \S~\ref{sec:coalesces} :
\begin{enumerate}
\item sample a random number of potential offspring using \S~\ref{sec:randomX}
\item record the surviving offspring each pointing to (or labelled
with) its immediate ancestor and add the labels to the tree
\S~\ref{sec:addgener} 
\item sample a new random sample \S~\ref{sec:randomsample}
\end{enumerate}
\item when a complete sample tree is found  the ancestry of the sample  is fixed so trace the
ancestry and record the branch lengths \S~\ref{sec:update},  \S~\ref{sec:recordonebls}
\item given the branch lengths for one sample  update the estimate of
$\EE{ \widetilde R_{i}^{N}(n)}$ \S~\ref{sec:ERiupdate}
\end{enumerate}



@*1 {\bf Includes}. 
\label{sec:includes}

The included libraries; we use the GSL library

@<includes@>=@#
#include <iostream>
#include <fstream>
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
#include <chrono>
#include <forward_list>
#include <assert.h>
#include <math.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>



@*1 {\bf Standard Library random number generator}.
\label{sec:stdlrng}

define the standard library random number generator

@<stdl rng@>=@#
/* \newline \aths{  obtain a seed out of thin air for the random number engine} */
  std::random_device randomseed;
  /* \newline \aths{  Standard Mersenne twister  random number engine}  */
  std::mt19937_64 rng(randomseed());


@*1 {\bf the parameters}. 
\label{sec:parameters}

define the parameters of the model, the   population size, the upper
bound $\zeta(N)$ and   $\alpha$  in \eqref{eq:2},  the
sample size $n$ and  the number of  experiments.     

@<parameters@>=@#
/* \newline \aths{ the population size $N$} */
const std::size_t CONST_POP_SIZE = 1e2; @#
/* \newline \aths{ the upper bound $\zeta(N)$ \eqref{eq:2}} */
const std::size_t CONST_CUTOFF = CONST_POP_SIZE; @#
/* \newline \aths{  $\alpha$ \eqref{eq:1}, \eqref{eq:2} }*/
const double CONST_ALPHA = 1.01; @#
/* \newline \aths{ sample size $n$} */
const std::size_t CONST_SAMPLE_SIZE = 1e1 ; @#
/* \newline \aths{ number of experiments}  */
const  int CONST_EXPERIMENTS = 1e2 ;


@*1 {\bf the gsl random number generator}.
\label{sec:gslrng}

define the GSL random number generator |rngtype| and the initialising function

@<gsl rng@>=@#
gsl_rng * rngtype ;
static void setup_rng( unsigned long int s )
{
        const gsl_rng_type *T ; 
        gsl_rng_env_setup(); 
        T = gsl_rng_default ;
        rngtype = gsl_rng_alloc(T);
        gsl_rng_set( rngtype,  s) ;
}



@*1 {\bf the probability mass function \eqref{eq:2} }.
\label{sec:pmf}

compute the kernel $k^{-\alpha} - (1+k)^{-\alpha}$ of  the    mass function in \eqref{eq:2}

@<pmf@>=@#
static double  kernel( const std::size_t &k )
{
  return ( pow(1./static_cast<double>(k), CONST_ALPHA) - pow( 1./static_cast<double>(k+1), CONST_ALPHA) ) ;
}


@*1 {\bf generate the CMF for sampling potential offspring }.
\label{sec:cmf}

generate the cumulative mass function (CMF) for sampling a random
number of potential offspring according to \eqref{eq:2}

@<cmf@>=@#
static void generatecmf( std::vector<double>& cmf)
{

@#

  double s {} ;
  for( std::size_t i = 1; i <= CONST_CUTOFF; ++i){
  /* \newline \aths{ |kernel| \S~\ref{sec:pmf}} */
  cmf[i] = cmf[i-1] + kernel( i);
    s += kernel(i); }

  assert(s > 0.);

/* \newline \aths{ normalise to generate a probability distribution} */

  std::transform( cmf.begin(), cmf.end(), cmf.begin(), [&s](const auto &x){return x/s;} );
  cmf[CONST_CUTOFF] = 1. ;
}


@*1 {\bf sample a random number of potential offspring}.
\label{sec:randomX}

sample a random number of potential offspring
$\min\{j : u \le F(j) \}$ where $u$ is a random uniform and $F$ the
cumulative mass function \S~\ref{sec:cmf}

@<sample a number of offspring@>=@#
static std::size_t randomX( const  std::vector<double>& f )
{

@#

  /* \newline \aths{ $f$ is the cumulative mass function generated in
  \S~\ref{sec:cmf}} */
  const double u = gsl_rng_uniform_pos( rngtype);
  std::size_t j {1} ;
 
  while( u > f[j]){ ++j ;}

  assert( j >= 1);
  return j ;
}

@*1 {\bf add a generation}.
\label{sec:addgener}

add a generation to the population tree by sampling $N$ surviving
offspring among the potential offspring produced by the current
individuals and recording the immediate ancestors of the new offspring

@<add a generation@>=@#
static void  addgeneration( std::vector<std::size_t>& tree, const std::vector<double>& vcmf )
{

 @#

  /* \newline \aths{  |y| records the pool of potential offspring}  */
  std::vector< std::size_t > y {} ;
  y.clear() ;
  std::size_t x {} ;
  
  for( std::size_t i = 0 ; i < CONST_POP_SIZE; ++i){
  /* \newline  \aths{ |randomX|  \S~\ref{sec:randomX}} */
    x = randomX(vcmf) ;
    y.reserve(y.size() +  x);

/* \newline \aths{  the individual on level $i$ produces $x$ potential
offspring so in the pool of potential offspring there will be $x$
offspring pointing to level $i$ } */

   y.insert( y.end(), x, i); }

   assert(y.size() >= CONST_POP_SIZE); 
  /* \newline  \aths{ add new generation to tree} */
  /* \newline \aths{ |y| is the immediate ancestors of the new offspring} */

/* \newline \aths{ shuffle the levels of the potential offspring and the
first $N$ will survive and be inserted into the tree, the population
ancestry } */

  std::shuffle( y.begin(), y.end(), rng); 
  tree.reserve( tree.size() + CONST_POP_SIZE) ;
  std::move( y.begin(), y.begin() + CONST_POP_SIZE, std::back_inserter(tree));

  y.clear();
    y.shrink_to_fit();
    std::vector<std::size_t>().swap(y);
}

@*1 {\bf  get a random sample}.
\label{sec:randomsample}

get a random sample by sampling labels uniformly at random without
replacement; a sample with $m$ blocks at time $g$ is the vector
$\left( (s_{1}, \ell_{1}), \ldots, (s_{m}, \ell_{m}) \right)$ where
$s_{i}$ is the size of block $i$ and $\ell_{i}$ the level of the
block; the immediate ancestor of the individual on level $\ell_i$ at
time $g$ is $A_{\ell_{i}}(g)$


@<random sample@>=@#
static void randomsample( std::vector< std::pair< std::size_t, std::size_t>>& sample )
{
  /* \newline  \aths{ pair is (size of block, level  of block)} */
  /* \newline \aths{ |V| will be the levels of the new sample} */
  std::vector< std::size_t> V (CONST_POP_SIZE, 0);
  std::iota( V.begin(), V.end(), 0);
  std::shuffle( V.begin(), V.end(), rng);
  sample.clear();
  assert( sample.size() < 1);

  sample.reserve(  CONST_SAMPLE_SIZE  ); 
  for( std::size_t i = 0 ; i < CONST_SAMPLE_SIZE ; ++i){
    /* \newline \aths{  pair is (size of block, level of block)} */
    /* \newline  \aths{ |V[i]| is the sampled  level of leaf |i|} */ 
    sample.push_back( std::make_pair( 1, V[i] ) ) ; }
  assert(sample.size() == CONST_SAMPLE_SIZE);
}


@*1 {\bf the immediate  ancestor}.
\label{sec:agi}

get the immediate ancestor (parent) $A_{\ell}(g)$ of the individual
living on level $\ell$ at time $g$

@<agi@>=@#
static std::size_t getagi( const std::size_t &g, const std::size_t &level,  const std::vector<std::size_t>& tree )
{

@#

  /* \newline  \aths{ get  $A_{|level|}(g)$;} */
  /* \newline \aths{  the immediate  ancestor of individual on level |level| at
  time  |g| } */

  return ( tree[ (g*CONST_POP_SIZE) + level]) ;  
}

@*1 {\bf check if the sample tree is complete}. 
\label{sec:coalesces}

check if the sample tree is complete with the leaves finding a common
ancestor; here the sampled leaves go searching  for a common ancestor

@<coalesced@>=@#
static bool allcoalesced( const std::size_t& generations, const std::vector<std::pair< std::size_t, std::size_t>>& sample, const std::vector< std::size_t>& tre )
{

/* \newline  |generations| + 1 is the current number of generations in
the tree (numbered from zero) */

  std::size_t g = generations; 

  std::vector< std::size_t> a( CONST_SAMPLE_SIZE, 0);
  /* \newline \aths{  pair is (size of block, level of block)} */
  /* \newline \aths{  copy sampled levels  into |a|; use |a| to check if the tree is complete}  */
   std::transform( sample.begin(), sample.end(), a.begin(), [](const auto &x){return x.second;});
  
  while( (a.size() > 1) && (g > 0)){
  /* \newline \aths{ record the immediate ancestors} */
    for( std::size_t i = 0 ; i < a.size() ; ++i){
      a[i] = getagi( g, a[i], tre); }
      /* \newline \aths{ remove duplicate entries signalling common ancestors} */
    std::sort( a.begin(), a.end());
    a.erase( std::unique( a.begin(), a.end()), a.end() );
    --g ;
  }
  /* \newline \aths{ return TRUE if sample tree is complete}  */
  return  (a.size() < 2) ;

}

@*1 {\bf compare size of block}.
\label{sec:compare}

compare size of block for  sorting in descending order on size of
block 

@<compare@>=@#
static bool comp( const std::pair< std::size_t, std::size_t> &a, const std::pair< std::size_t, std::size_t> &b )
{
  return a.first > b.first ;
}

@*1 {\bf update sample}.
\label{sec:update}

update a  sample given that the tree is complete; merge blocks and record
continuing  blocks 

@<update@>=@#
static void updatesample( const std::size_t& g, std::vector<std::pair< std::size_t, std::size_t>>& sample, const std::vector<std::size_t>& tree )
{
  /* \newline  \aths{  pair is (size of block, level of  block)} */
  
   std::size_t s {} ;
   const std::size_t e = sample.size();
   std::size_t k {} ;

   for( std::size_t i = 0 ; i < e; ++i){
    s = 0 ;
    for( std::size_t j = i ; j < e ; ++j){
      s += ( getagi(g, sample[i].second ,tree) == getagi(g, sample[j].second, tree) ? sample[j].first : 0) ;
      /* \newline  \aths{ if block has already been merged set size of block to zero} */
      sample[j].first = ( getagi(g,  sample[i].second, tree) == getagi(g, sample[j].second, tree) ? 0 : sample[j].first); } 
    if( s > 0 ){
      /* \newline \aths{  record only a new merged block or a current  block not merging} */
      ++k ;
      /* \newline \aths{  record the continuing block with the block size} */
      /* \newline  \aths{ and the level of the immediate ancestor of the block} */
      sample.push_back( std::make_pair( s,  getagi( g, sample[i].second, tree) ) ) ; }}

/* \newline  \aths{  sort the sample on block size in descending order using
\S~\ref{sec:compare}}  */
   std::sort( sample.begin(), sample.end(), comp);
   /* \newline \aths{ remove all blocks with block size zero}  */
   sample.resize(k);
   assert( sample.back().first > 0);

   assert( sample.size() > 1 ? (sample.back().first < CONST_SAMPLE_SIZE) :   (sample.back().first == CONST_SAMPLE_SIZE) ); 
}


@*1 {\bf update the current spectrum}. 
\label{sec:updatecurrentbls}

update the current branch length spectrum

@<current bls@>=@#
static void updatevbi( const std::vector< std::pair< std::size_t, std::size_t>>& sample, std::vector<std::size_t>& vbi )
{
  /* \newline  \aths{  pair is (size of block, level of block)} */
  for( std::size_t i = 0; i < sample.size(); ++i){
  /* \newline \aths{ update the total tree length} */
    ++vbi[0] ;
    assert( sample[i].first > 0);
    assert( sample[i].first < CONST_SAMPLE_SIZE );
    /* \newline  \aths{  update  $\ell_{j}^{N}(n)$ where |j = sample[i].first|} */
    ++vbi[ sample[i].first] ; }
}


@*1 {\bf record one branch length spectrum}.
\label{sec:recordonebls}

record one realisation of
$\left (L_{1}^{N}(n), \ldots, L_{n-1}^{N}(n) \right)$; the sample tree
is complete so read the branch lengths off the tree



@<one bls@>=@#
static void recordonebls( std::vector< std::pair< std::size_t, std::size_t>>& csample,  const std::vector<std::size_t>& tre, std::vector<std::size_t>& bls )
{

  std::size_t g =  (tre.size()/CONST_POP_SIZE)  - 1   ;
  std::fill( bls.begin(), bls.end(),0);
  while( csample.back().first < CONST_SAMPLE_SIZE ){
  /* \newline  \aths{ \S~\ref{sec:updatecurrentbls} for |updatevbi| } */
    updatevbi( csample, bls);
    /* \newline \aths{  \S~\ref{sec:update} for |updatesample|} */
    updatesample( g, csample, tre);
    --g ;}
}


@*1 {\bf update the estimate $\overline\rho_{i}^{N}(n)$ of $\EE{\widetilde R_{i}^{N}(n)}$}.
\label{sec:ERiupdate}

@<$\overline\rho_{i}^{N}(n)$ update@>=@#
static void updaterestimate( const std::vector< std::size_t>& b, std::vector<double>& r)
{
/* \newline \aths{ |b[0]| is total tree size} */
  assert( b[0] >=  CONST_SAMPLE_SIZE );
  for( std::size_t i = 1; i < CONST_SAMPLE_SIZE; ++i){
    r[i] += static_cast<double>( b[i])/static_cast<double>(b[0]) ; }
}


@*1 {\bf one ancestry}. 
\label{sec:ancestry}

record one ancestry and the resulting branch lengths; add to the
population tree and draw samples until a sample tree is complete, then
trace the ancestry of the sample and record the branch lengths


@<ancestry@>=@#
static void oneancestry( std::vector<double>& vr, const std::vector<double>& vf)
{

@#

  std::vector< std::size_t> tre (CONST_POP_SIZE, 0);
  std::iota( tre.begin(), tre.end(), 0);
  std::vector< std::pair< std::size_t, std::size_t>> v {} ;
  std::vector< std::size_t> b (CONST_SAMPLE_SIZE, 0);

  /* \newline \aths{  \S~\ref{sec:randomsample} for |randomsample|} */
  randomsample(v);
  std::size_t numbergenerations = 0 ;
  /* \newline  \aths{  \S~\ref{sec:coalesces} for |allcoalesced|} */
  while( !allcoalesced( numbergenerations, v, tre)){
  /* \newline \aths{ \S~\ref{sec:addgener}} */
    addgeneration( tre, vf);
    randomsample(v);
    ++numbergenerations; }

    /* \newline  \aths{ sample merges so record the resulting bls;
    \S~\ref{sec:recordonebls} for |recordonebls| } */
  recordonebls( v, tre, b);

  /* \newline \aths{  update estimate of relative branch lengths;
  \S~\ref{sec:ERiupdate} for |updaterestimate|} */
  updaterestimate( b, vr) ;
}


@*1 {\bf approximate  $\EE{\widetilde R_{i}^{N}(n)}$}.
\label{sec:estimate}


estimate  $\EE{\widetilde R_{i}^{N}(n)}$ by generating |CONST_EXPERIMENTS|
number of ancestries and tracing each for the branch lengths; recall
\S~\ref{sec:parameters} for the parameter values

@<estimate@>=@#
static void estimate()
{
  /* \newline \aths{ estimate the relative branch lengths} */
  std::vector<double> vcmfx (CONST_CUTOFF + 1, 0 ) ;
  std::vector< double> vR (CONST_SAMPLE_SIZE, 0);

  /* \newline \aths{ \S~\ref{sec:cmf} for |generatecmf|} */
  generatecmf( vcmfx); 
  int r = CONST_EXPERIMENTS + 1 ;
  while( --r > 0){
  /* \newline \aths{ \S~\ref{sec:ancestry} for |oneancestry|} */
    oneancestry(vR, vcmfx); }

/* \newline \aths{ record the estimate  $\overline\rho_{i}^{N}(n)$   of $\EE{\widetilde R_{i}^{N}(n)}$} */
   for( const auto & z: vR){ std::cout << z/static_cast<double>(CONST_EXPERIMENTS) << '\n';}
}




@*1 {\bf the main module}. 
\label{sec:main}


The |main| function 

@C

/* \newline \S~\ref{sec:includes} */
@<includes@>@#
/* \newline \S~\ref{sec:stdlrng} */
@<stdl rng@>@#
/* \newline \S~\ref{sec:parameters} */
@<parameters@>@#
/* \newline \S~\ref{sec:gslrng} */
@<gsl rng@>@#
/* \newline \S~\ref{sec:pmf} */
@<pmf@>@#
/* \newline \S~\ref{sec:cmf} */
@<cmf@>@#
/* \newline \S~\ref{sec:randomX} */
@<sample a number of offspring@>@#
/* \newline \S~\ref{sec:addgener} */
@<add a generation@>@#
/* \newline \S~\ref{sec:randomsample} */
@<random sample@>@#
/* \newline \S~\ref{sec:agi} */
@<agi@>@#
/* \newline \S~\ref{sec:coalesces} */
@<coalesced@>@#
/* \newline \S~\ref{sec:compare} */
@<compare@>@#
/* \newline \S~\ref{sec:update} */
@<update@>@#
/* \newline \S~\ref{sec:updatecurrentbls} */
@<current bls@>@#
/* \newline \S~\ref{sec:recordonebls} */
@<one bls@>@#
/* \newline \S~\ref{sec:ERiupdate} */
@<$\overline\rho_{i}^{N}(n)$ update@>@#
/* \newline \S~\ref{sec:ancestry} */
@<ancestry@>@#
/* \newline \S~\ref{sec:estimate} */
@<estimate@>@#


int main(int argc, const  char * argv[])
{
/* \newline \aths{ \S~\ref{sec:gslrng} for |setup_rng|} */
  setup_rng( static_cast<unsigned long>( atoi(argv[1]))) ;
/* \newline \aths{ \S~\ref{sec:estimate} for |estimate|} */
  estimate() ;
  gsl_rng_free(rngtype);
  return GSL_SUCCESS ;
}



@* {\bf examples}.
\label{sec:examples}

Let the quenched $(N,n)$-coalescent denote the trees generated by
evolving a haploid population forward in time until a random sample of
$n$ leaves coalesces (finds a most recent common ancestor); the sample
then has a fixed tree or ancestry (recall \S~\ref{sec:intro}).  Since
we do not at this time have a quenched coalescent as $N\to \infty$ we
compare the normalised branch lengths $\EE{\widetilde R_{i}^{N}(n)}$
predicted by the quenched $(N,n)$-coalescent to $\EE{R_{i}(n)}$
predicted by the annealed coalescent obtained in the usual way as
$N\to \infty$.



In Fig~\ref{fig:quenchedannealedsfs} we compare the branch length
spectrum  $\EE{\widetilde R_{i}^{N}(n)}$ predicted by the
quenched $(N,n)$-coalescent  to  $\EE{R_{i}(n)}$  predicted by the annealed
incomplete $n$-Beta-coalescent  with coalescent rates when $k =
2,3,\ldots, n$ 
\begin{equation}
\label{eq:3}
\lambda_{n,k} =  \frac{1}{B(\gamma;  2-\alpha,\alpha)} \int_{0}^{1}\one{0 < x \le \gamma } x^{k-\alpha - 1}(1-x)^{n+\alpha - k - 1} dx 
\end{equation}
where $B(\gamma ; 2-\alpha,\alpha) =
\int_{0}^{1}\one{0 < x\le \gamma }x^{1-\alpha}(1-x)^{\alpha-1}dx $ and
\begin{equation}
\label{eq:4}
\gamma =  \frac{K}{K + {m}_{\infty}}
\end{equation}
when $\zeta(N) = KN$ recall \eqref{eq:2}  for some constant $K > 0$, 
and $m_{\infty} =  \lim_{N\to \infty}\EE{X_{1}} $ the
expected number of potential offspring in an arbitrarily large
population; we approximate  ${m}_{\infty}$    with 
\begin{equation}
\label{eq:5}
  {m}_{\infty} \approx  \frac 12 \left(2 +  \frac{ 1 + 2^{1-\alpha}}{\alpha - 1} \right)
\end{equation}


\begin{figure}[htp]
\centering
\subfloat[$\zeta(N) = \infty$]{\includegraphics[scale=0.4]{quenchedunboundedannealedbetasfsNe3-crop}}
\subfloat[$\zeta(N) = N$]{\includegraphics[scale=0.4]{quenchedannealedsfsNe3-crop}}\\
\caption[quenched vs.\ annealed]{Comparing
$\EE{\widetilde R_{i}^{N}(n)}$ and $\EE{R_{i}(n) }$.
Relative branch lengths compared between the quenched coalescent
(symbols) for $N=10^{3}$, $\alpha = 1.05$, $\zeta(N) = \infty$
according to \eqref{eq:1} (a) $\zeta(N) = N$ according to
\eqref{eq:2} (b) compared to the BLS predicted by the annealed
complete (a) and incomplete (b) Beta-coalescent \eqref{eq:3} (red
lines) for sample size $n$ as shown with $\gamma$ as in
\eqref{eq:4} where ${m}_{\infty}$ as in \eqref{eq:5} , recall we
take $\zeta(N) = N$ for the quenched coalescent. Here we graph
$\log(r_{i}(n)) - \log(1 - r_{i}(n))$ as a function of
$\log(i/n) - \log(1 - i/n)$ for $i = 1,2, \ldots, n-1$ where
$r_{i}(n)$ is an estimate of
$\EE{\widetilde R_{i}^{N}(n)}$ (symbols) and
$\EE{R_{i}(n) }$ (red lines) recall the notation from
\S~\ref{sec:intro}; the estimates of
$\EE{\widetilde R_{i}^{N}(n)}$  resp.\ $\EE{R_{i}(n) }$
from $10^{4}$ resp.\ $10^{6}$ experiments }
\label{fig:quenchedannealedsfs}
\end{figure}


@* {\bf conclusions and bibliography}. 
\label{sec:concl}

Write $(x)_{m} = x(x-1)\cdots (x-m+1)$ for any real $x$ and $m\in \mathds N$ and $(x)_{0} \equiv 1$.    
Let $\nu_{1},\ldots, \nu_{N}$ denote the random number of offspring of
individuals in a haploid panmictic population of constant size $N$  current in an arbitrary generation, $\sum_{i}\nu_{i} = N$.
Let $c_{N} \equiv \EE{\nu_{1}(\nu_{1}-1) }/(N-1)$. If $c_{N}\to 0$ and 
\begin{equation}
\label{eq:limits}
\lim_{N\to \infty} \frac{\EE{ (\nu_{1})_{k_{1}} \cdots   (\nu_{r})_{k_{r}}  } }{ N^{k_{1} + \cdots + k_{r} - r}  c_{N} }
\end{equation}
exist for all $k_{1}, \ldots, k_{r} \ge 2 $ and $r \in \mathds N$ then
$\set{\xi^{n,N}\svigi{\lfloor t/c_{N} \rfloor } ; t \ge 0 }$ converge
to $\set{\xi^{n}(t); t\ge 0 }$ (in finite-dimensional distributions)
and the transition rates of $\set{\xi^{n}(t); t\ge 0 }$ are uniquely
determined by the limit in \eqref{eq:limits}
\citep[Proposition~1]{schweinsberg03} and \cite{MS01}.  However,
\eqref{eq:limits} implicitly averages over the ancestral relations of
the gene copies in a sample.

With the C++ code presented here one can use simulations to check if
conditioning on the (random) population ancestry   matters for
predictions of genetic variation, and thus for inference of 
evolutionary histories of natural populations. 


%%\bibliographystyle{alpha}
%%\bibliography{refs}

\begin{thebibliography}{DFBW24}

\bibitem[DFBW24]{Diamantidis2024}
Dimitrios Diamantidis, Wai-Tong~(Louis) Fan, Matthias Birkner, and John
  Wakeley.
\newblock Bursts of coalescence within population pedigrees whenever big
  families occur.
\newblock {\em GENETICS}, 227(1), February 2024.

\bibitem[KL94]{knuth1994cweb}
Donald~Ervin Knuth and Silvio Levy.
\newblock {\em The CWEB system of structured documentation: version 3.0}.
\newblock Addison-Wesley Longman Publishing Co., Inc., Reading, Massachusetts,
  1994.

\bibitem[KR88]{kernighan1988c}
Brian~W Kernighan and Dennis~M Ritchie.
\newblock The {C} programming language, 1988.

\bibitem[MS01]{MS01}
M~M\"{o}hle and S~Sagitov.
\newblock A classification of coalescent processes for haploid exchangeable
  population models.
\newblock {\em Ann Probab}, 29:1547--1562, 2001.

\bibitem[Sch03]{schweinsberg03}
J~Schweinsberg.
\newblock Coalescent processes obtained from supercritical {G}alton-{W}atson
  processes.
\newblock {\em Stoch Proc Appl}, 106:107--139, 2003.

\bibitem[Tan11]{tange11:_gnu_paral}
O~Tange.
\newblock {GNU} parallel -- the command-line power tool.
\newblock The USENIX Magazine, 2011.

\end{thebibliography}





@
\end{document}
