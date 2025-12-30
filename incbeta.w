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
%%\usepackage[round,numbers,super]{natbib}
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
\setstretch{1}
%\usepackage{abstract}
%\renewcommand{\abstractname}{}
\newcommand{\one}[1]{\ensuremath{\mathds{1}_{\left\{ #1 \right\}}}}
\newcommand{\EE}[1]{\ensuremath{\mathds{E}\left[ #1 \right]}}%
\newcommand{\im}{\ensuremath{\imath} }%
\newcommand{\jm}{\ensuremath{\jmath} }%
\newcommand{\be}{\begin{equation}}%
\newcommand{\ee}{\end{equation}}%
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
  {\textsf{\textbf{\@@title}}}
\end{flushleft}

\begin{center}
\textsc{\@@author}
\end{center}
\egroup
}
\makeatother
\title{Beta coalescents when sample size is large \\ --- estimating $\EE{R_{i}(n)}$ for  the  Beta$(\gamma,2-\alpha,\alpha)$-coalescent}
\author{Bjarki Eldon\footnote{\href{beldon11@@gmail.com}{beldon11@@gmail.com}} \footnote{Supported by 
  Deutsche Forschungsgemeinschaft (DFG) - Projektnummer 273887127 
%% „funded by the Deutsche Forschungsgemein-schaft (DFG, German Research Foundation) –Projektnummer(n)“.
%% 
through DFG SPP 1819: Rapid Evolutionary Adaptation grant STE 325/17
to Wolfgang Stephan; acknowledge funding by the Icelandic Centre of
Research (Rann\'is)  through an Icelandic Research Fund (Ranns\'oknasj\'o{\dh}ur)  Grant of Excellence no.\
185151-051 to Einar \'Arnason, Katr\'in Halld\'orsd\'ottir, Alison M.\
Etheridge, Wolfgang Stephan, and BE;  Start-up
module grants through SPP 1819 with Jere Koskela and Maite
Wilke-Berenguer, and with Iulia Dahmer. \\ \today }
\orcidlink{0000-0001-9354-2391} }

\begin{document}
\maketitle
\renewcommand{\abstractname}{\vspace{-\baselineskip}}


%%\rule{\textwidth}{.8pt}


\begin{abstract}
  Let  $\{ \xi^{n}(t) : t \ge 0 \}$ be 
the Beta$(\gamma,2-\alpha,\alpha)$-coalescent, $\# A$ is
the cardinality of a given set $A$, $n$ sample size,  
$L_{i}^{N}(n) \equiv \int_{0}^{\tau(n) } \# \left\{ \xi \in \xi^{n}(t)
: \#\xi = i \right\}dt $ and   $L(n) \equiv \int_{0}^{\tau(n)} \# \xi^{n}(t)dt $  and 
$\tau(n) \equiv \inf \left\{ t \ge 0 : \# \xi^{n}(t) = 1 \right\} $
for $i \in \{1, 2, \ldots, n-1\}$;
$R_{i}(n) \equiv L_{i}(n)/\sum_{j}L_{j}(n) $ for $i=1,2,\ldots, n-1$.
Then $L_{i}^{N}(n)$ is interpreted as the random total length of
branches supporting $i \in \{1, 2, \ldots, n-1\}$ leaves, with the
length measured in coalescent time units,  and $n$ sample size.  We then have
$L(n) = L_{1}(n) + \cdots + L_{n-1}(n)$.
With this C++ code one estimates the functionals 
$\EE{R_{i}(n)}$ of gene genealogies described by the
Beta$(\gamma,2-\alpha,\alpha)$-coalescent where $0 < \gamma \le 1$ and
$1 \le \alpha < 2$. The Beta$(\gamma,2-\alpha,\alpha)$-coalescents are
a family of $\Lambda$-coalescents \cite{P99,DK99,S99}; the transition
rates are
\begin{displaymath}
\lambda_{n,k} =   \binom{n}{k}\frac{ B(\gamma,k-\alpha,n-k+\alpha)}{B(\gamma,2-\alpha,\alpha)}
\end{displaymath}
for $k = 2, 3, \ldots, n$ where $B(x,a,b) = \int_{0}^{x}
t^{a-1}(1-t)^{b-1}dt$ for $0 < x \le 1$ and $a,b > 0$
\cite{chetwyn-diggle_beta}.  The
Beta$(\gamma,2-\alpha,\alpha)$-coalescent extends the
Beta$(2-\alpha,\alpha)$-coalescent derived from a model of
sweepstakes reproduction  (skewed offspring number distribution)
\cite{schweinsberg03}. 
\end{abstract}

\tableofcontents


@* {\bf Copyright}. 


Copyright {\textcopyright} {\the\year}  Bjarki Eldon \newline


{\tt incbeta : simulate branch lengths from an incomplete beta-coalescent}




This document and any source code it contains is distributed under the
terms of the GNU General Public License (version $\ge 3$).  You should
have received a copy of the license along with this file (see file
COPYING).

    The source codes described in this document are free software: you
    can redistribute it and/or modify it under the terms of the GNU
    General Public License as published by the Free Software
    Foundation, either version 3 of the License, or (at your option)
    any later version.

    This document and the code it contains   is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this file (see COPYING).  If not, see \url{http://www.gnu.org/licenses/}.


@* {\bf Compilation,  output and execution}. 
\label{compile}

 This CWEB
      \cite{knuth1994cweb} document (the {\tt .w} file) can be
      compiled with {\tt cweave} to generate a {\tt .tex} file, and
      with {\tt ctangle} to generate a {\tt .c} \cite{kernighan1988c}
      file.

One can use {\tt cweave} to generate a {\tt .tex} file, and {\tt
ctangle} to generate a {\tt .c} file. To compile the C++ code (the {\tt
.c} file), one needs the GNU Scientific Library.


Compiles on {\tt Linux Debian trixie/sid}  with kernel {\tt 6.12.10-amd64}  and   {\tt ctangle 4.11} and {\tt g++ 14.2} and {\tt GSL 2.8}


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

 Let  $\{ \xi^{n}(t) : t \ge 0 \}$ be 
the Beta$(\gamma,2-\alpha,\alpha)$-coalescent, $\# A$ is
the cardinality of a given set $A$, $n$ sample size,  
$L_{i}^{N}(n) \equiv \int_{0}^{\tau(n) } \# \left\{ \xi \in \xi^{n}(t)
: \#\xi = i \right\}dt $ and $L(n) \equiv \int_{0}^{\tau(n)} \# \xi^{n}(t)dt $ and 
$\tau(n) \equiv \inf \left\{ t \ge 0 : \# \xi^{n}(t) = 1 \right\} $
for $i \in \{1, 2, \ldots, n-1\}$;
$R_{i}(n) \equiv L_{i}(n)/\sum_{j}L_{j}(n) $ for $i=1,2,\ldots, n-1$.
Then $L_{i}^{N}(n)$ is interpreted as the random total length of
branches supporting $i \in \{1, 2, \ldots, n-1\}$ leaves, with the
length measured in coalescent time units,  and $n$ sample size.  We then have
$L(n) = L_{1}(n) + \cdots + L_{n-1}(n)$.


  We estimate   $\EE{R_{i}(n)}$ when 
the gene genealogy is determined by the
Beta$(\gamma,2-\alpha,\alpha)$-coalescent with transition rates
\begin{equation}
\label{eq:lambdank}
\lambda_{n,k} =  \binom{n}{k}\frac{ B(\gamma,k-\alpha,n-k+\alpha)}{B(\gamma,2-\alpha,\alpha)}
\end{equation}
for $2 \le k \le n$ and
$B(\gamma,a,b) = \int_{0}^{\gamma} t^{a-1}(1-t)^{b-1}dt$ and $a,b > 0$
and $0 < \gamma \le 1$.

The Beta$(\gamma,2-\alpha,\alpha)$-coalescent can be shown to describe
the random gene genealogies of a sample when the sample comes from a
haploid panmictic population of constant size evolving according to 
randomly increased recruitment and when there is an upper bound on the random 
number of potential offspring any arbitrary individual can produce.
The upper bound translates to the parameter $\gamma$ of the
Beta$(\gamma,2-\alpha,\alpha)$-coalescent \cite{chetwyn-diggle_beta}.


The algorithm is summarised in \S~\ref{sec:code}, the code follows 
in \S~\ref{sec:includes}--\S~\ref{sec:main}; we conclude in
\S~\ref{sec:concl}. Comments within the code are in \aths{this font and colour}



@* {\bf Code}. 
\label{sec:code}


Write $[n] = \{1,2, \ldots, n\}$ for any natural number $n$.  Let $n$
denote the sample size and $m \in [n]$ the current number of blocks.
We are  interested in the branch lengths and  require the
block sizes $(b_{1}, \ldots, b_{m})$ where $b_{j} \in [n]$ and
$b_{1} + \cdots + b_{m} = n$ are the current block  sizes

\begin{enumerate}
\item $\left( r_{1}(n), \ldots, r_{n-1}(n)\right) \leftarrow (0,\ldots, 0)$
\item for each of  $M$ experiments; \S~\ref{sec:estimateri}:
\begin{enumerate}
\item set $(b_{1}, \ldots, b_{n}) \leftarrow  (1,\ldots, 1)$.  
\item set current branch lengths $\ell_{i}(n) \leftarrow  0$ for $i \in [n-1]$
\item set the current number of blocks  $m\leftarrow n$ 
\item {\bf while} $m > 1$ (at least two blocks; see \S~\ref{sec:onetree})  :
\begin{enumerate}
\item sample exponential time $t$ with rate  $\lambda_{m, 2} + \cdots + \lambda_{m,m}$ \eqref{eq:lambdank}
\item update the current  branch lengths   $\ell_{b}(n) \leftarrow  t +  \ell_{b}(n)$ for $b = b_{1},\ldots, b_{m}$; \S~\ref{sec:updateLin}
\item sample merger size $j \in \{2,3,\ldots, m\}$  \S~\ref{sec:mergersize}
\item given merger size shuffle the blocks and merge $j$ of  them   \S~\ref{sec:onetree}
\item update current number of blocks  $m \leftarrow m - j + 1$
\end{enumerate}
\item given a realisation $\ell_{1}(n), \ldots, \ell_{n-1}(n)$   of branch lengths update the estimate of
$\EE{R_{i}(n)}$; $r_{i}(n) \leftarrow r_{i}(n) + \ell_{i}(n)/\sum_{j}\ell_{j}(n)$  \S~\ref{sec:updateri}
\end{enumerate}
\item return an estimate  $(1/M)r_{i}(n)$  of   $\EE{R_{i}(n)}$ for $i = 1,2, \ldots, n-1$
\end{enumerate}




@*1 {\bf includes}. 
\label{sec:includes}

the included libraries; we use the {\tt GSL} and {\tt  boost}   libraries

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
#include <gsl/gsl_sf.h>
#include <boost/math/special_functions/beta.hpp>


@*1 {\bf the random number generator}. 
\label{SEC:rng}


initialise the random number engines; we will not go into discussions
about how to get a computer to give us a ``random'' number


@<gslrng@>=@#
/* \newline  \aths{the GSL random number engine }  */
gsl_rng * rngtype ;
 /* \newline \aths{  obtain a seed out of thin air for the STL  random number engine }  */
 std::random_device randomseed;
  /* \newline  \aths{  The STL  standard Mersenne twister   random number engine seeded with |randomseed()| }
    */
  std::mt19937_64 rng(randomseed());

/* \newline  \aths{ set up and initialise the GSL random number generator }  */
static void setup_rng( unsigned long int s )
{
        const gsl_rng_type *T ; 
        gsl_rng_env_setup(); 
        T = gsl_rng_default ;
        rngtype = gsl_rng_alloc(T);
        gsl_rng_set( rngtype,  s) ;
}


@*1 {\bf incomplete beta function}.
\label{SEC:incbetaHGF}

use the GSL Gauss hypergeometric function to compute the incomplete
beta function; we have the representations
\begin{displaymath}
\begin{split}
B(x,a,b) &  =   x^{a} F(a, 1-b; a+1; x)/a  \\
B(x,a,b) &  =   x^{a}(1-x)^{b} F(a + b, 1; a+1; x) / a
\end{split}
\end{displaymath}

return the logarithm of $B(x,a,b)$ 

@<log of incomplete Beta@>=@#
static long double  lnincbetaGF(const long double& a, const long double &b, const long double &x)
{
   
  return  logl( static_cast<long double>( gsl_sf_hyperg_2F1( a+b, 1, a+1, x)) ) +  (a*logl( x )) + (b*log(1.-x)) - log(a) ;
}


@*1 {\bf the incomplete beta function using boost}.
\label{sec:incbetaboost}


the incomplete beta function using the boost library

@<incomplete beta using boost@>=@#
static double incbeta(  const   double &a, const  double &b,  const  double &x )
{
   assert( x <= 1.);
  assert( 0 <= x);
  /* \newline  \aths{ if using the GSL library    |gsl_sf_beta_inc( a, b, x) * gsl_sf_beta(  a, b  )| }  */
  
  
   return ( x < 1 ? boost::math::beta( static_cast<long double>(a), static_cast<long double>(b), static_cast<long double>(x) ) : gsl_sf_beta( a,b) )   ;
}


@*1 {\bf the merger rate}.
\label{sec:lambdank}


compute the merger rate $\binom{n}{k}B(\gamma,k-\alpha, n-k +\alpha)/B(\gamma,2-\alpha,\alpha)$

@<merger rate@>=@#
static double rate(  const   double &m, const   double &k, const  double &a,  const  long double &x )
{

@#

  assert( k-a > 0);
  assert( k <= m);
  assert( m+a -k > 0);
  /* \newline \aths{ using |incbeta| from   \S~\ref{sec:incbetaboost} } */
  return static_cast<double>( expl( lgammal( m + 1) - lgammal(k+1) - lgammal(m-k+1) + logl( incbeta(k-a, m+a-k, x)) - logl(incbeta(2-a,a,x)) ) ) ;      
}


@*1 {\bf the total merger rate}.
\label{sec:lambdan}

compute the total merger rate $\lambda_{n} = \lambda_{n,2} + \cdots + \lambda_{n,n}$


@<lambdan@>=@#
static void totalrate( const double &n, const double &a, const double &x,  std::vector<double>& v )
{

 @#

  for( double m = 2 ; m <= n ; ++m){
    for ( double j = 2 ; j <= m ; ++j){
      assert( j <= m);
      /* \newline \aths{ using |rate| from \S~\ref{sec:lambdank}} */
      v[m] += rate( m, j, a, x); }}
}


@*1 {\bf sample merger size}.
\label{sec:mergersize}


sample merger size using the transition rates, returning
$\min\{ j : U \le \sum_{i=2}^{j}\lambda_{n,i}/\lambda_{n} \}$ where
$U$ is a random uniform from the unit interval and
$\lambda_{n} \equiv  \lambda_{n,2} + \cdots + \lambda_{n,n}$ \eqref{eq:lambdank}


@<mergersize@>=@#
static double getmerger( const double &m, const double &a, const double &x, const std::vector<double>& v )
{

  /* \newline \aths{ |m| is the current number of lines; |a| is $\alpha$; |x| is $\gamma$; } */
  /* \newline \aths{ |v| stores the  $\lambda_{n}$ values} */



  /* \newline \aths{ sample a random uniform} */
  const double u = gsl_rng_uniform( rngtype);
  double j = 2;
  /* \newline \aths{ |rate| from \S~\ref{sec:lambdank}} */
  double s =  rate( m, j, a, x) ;

@#

  while( u > s/v[ static_cast<int>(m) ]){  ++j; assert( j <= m); s += rate( m, j, a, x);  }

  return j;
}


@*1 {\bf update branch lengths $L_{i}(n)$}.
\label{sec:updateLin}

update the branch lengths $L_{i}(n)$ from a current configuration of
block sizes 


@<updateLin@>=@#
static void updateb( const double & timi, const std::vector<int>& tre,  std::vector<double>& b )
{
 /* \newline \aths{ |timi| is the sampled waiting time in the configuration given in |tre| } */
  assert( timi > 0);

@#

  std::for_each( tre.begin(), tre.end(), [&timi, &b]( const int t ){ assert(t > 0);  b[0] += timi; b[t] += timi; } );
}



@*1 {\bf update estimate of $\EE{R_{i}(n)}$}.
\label{sec:updateri}

update the estimate of $\EE{R_{i}(n)}$ for $i = 1, 2, \ldots, n-1$


@<updateri@>=@#
static void updateri( const std::vector<double>& bi, std::vector<double>& ri)
{

@#

  const double d = bi[0] ;

@#
  assert( d > 0);

@#
  std::transform( bi.begin(), bi.end(), ri.begin(), ri.begin(), [&d]( const auto &x, const auto &y ){return y + (static_cast<double>(x)/d); });
}



@*1 {\bf one tree}.
\label{sec:onetree}


generate one realisation of $L_{i}(n)$ for  $i = 1, 2, \ldots, n-1$

@<generate one tree@>=@#
static void genealogy( const int &n, const double &a, const double &x,   const std::vector<double>& v, std::vector<double>& vri )
{

@#

  std::vector<int> t (n, 1);
  std::size_t ms {} ;
  double timi {} ;
  int newb {} ; 
  std::vector<double> vb( n) ;
  std::size_t q {} ;
  
  while(t.size() > 1){
  /* \newline \aths{  sample waiting  time until next merger}  */
    timi = gsl_ran_exponential( rngtype, 1./v[ t.size() ]);
    assert(timi > 0);
    /* \newline \aths{ update branch lengths \S~\ref{sec:updateLin}} */
    updateb( timi, t, vb);
    /* \newline \aths{  get the size of next merger \S~\ref{sec:mergersize} } */
    ms = static_cast<std::size_t>( getmerger( static_cast<double>( t.size() ), a, x, v) );
    /* \newline \aths{  shuffle the blocks and merge the rightmost |ms| blocks } */
    std::shuffle( t.begin(), t.end(), rng);
    /* \newline \aths{ get the size of the new block } */
    newb = std::accumulate( t.rbegin(), t.rbegin() + ms, 0 );
    /* \newline \aths{ |q| is the current number of blocks } */
    q = t.size() ;
    /* \newline \aths{ remove the merged blocks} */
    t.resize( q - ms );
    /* \newline \aths{ add the new block |newb|  to the configuration } */
    t.push_back( newb); }

/* \newline \aths{ given realised branch lengths  update the estimate of $\EE{R_{i}(n)}$ \S~\ref{sec:updateri}} */
  updateri( vb, vri);
}



@*1 {\bf estimate $\EE{R_{i}(n)}$}.
\label{sec:estimateri}


@<get an estimate of $\EE{R_{i}(n)}$@>=@#
static void estimate( const double &n, const double &a, const double& K)
{




/* \newline  \aths{ approximate  the mean $|mu| =  m_{\infty} \approx (2 + (1 + 2^{1-\alpha})/(\alpha - 1))/2$} */

@#

/* \newline  \aths{  need $\alpha > 1$ when applying a cutoff; see the approximation of  $m_{\infty}$ below } */ 
 const double mu = (a > 1 ? ((1 + (pow(2., 1. - a)/(a - 1))) + (1 + (1/(a - 1))) )/2. : 0) ;

@#

/* \newline \aths{  $K$ is the cutoff constant;  $K = 0$ is  taken
as unbounded distribution of number of potential offspring translating
to complete Beta$(2-\alpha,\alpha)$-coalescent; otherwise the cutoff is $K/(m_{\infty} + K)$} */

@#

/* \newline \aths{ if $\alpha = 1$ taking the cutoff as $K$} */
const double p = (a > 1 ?  (K > 0 ? K/( mu + K ) : 1) : K) ;

@#
  

  std::vector<double> v ( static_cast<int>(n) + 1 );
  /* \newline \aths{ |totalrate|  \S~\ref{sec:lambdan}} */
  totalrate( n, a, p, v);

  std::vector<double> vri (static_cast<int>(n) );
  /* \newline \aths{ set to $10^{5}$ number of experiments} */
  int r = 1e5 + 1;
  while( --r > 0){
  /* \newline \aths{  |genealogy|  \S~\ref{sec:onetree}} */
    genealogy( static_cast<int>(n), a, p, v, vri); }

/* \newline \aths{ print the estimates of $\EE{R_{i}(n)}$ summer over the experiments} */
  std::for_each( vri.begin(), vri.end(), [](const auto &x){ std::cout << x << '\n'; });
}





@*1 {\bf the main module}. 
\label{sec:main}

The |main| function

@C

/* \newline \S~\ref{sec:includes} */
@<includes@>@#
/* \newline \S~\ref{SEC:rng} */
@<gslrng@>@#
/* \newline \S~\ref{SEC:incbetaHGF} */
@<log of incomplete Beta@>@#
/* \newline \S~\ref{sec:incbetaboost} */
@<incomplete beta using boost@>@#
/* \newline \S~\ref{sec:lambdank} */
@<merger rate@>@#
/* \newline \S~\ref{sec:lambdan} */
@<lambdan@>@#
/* \newline \S~\ref{sec:mergersize} */
@<mergersize@>@#
/* \newline \S~\ref{sec:updateLin} */
@<updateLin@>@#
/* \newline \S~\ref{sec:updateri} */
@<updateri@>@#
/* \newline \S~\ref{sec:onetree} */
@<generate one tree@>@#
/* \newline \S~\ref{sec:estimateri} */
@<get an estimate of  $\EE{R_{i}(n)}$@>@#



int main(int argc, char *argv[])
{

/* \newline \aths{  initialise the GSL random number generator |rngtype|
 using  |setup_rng|  \S~\ref{SEC:rng} } */
setup_rng(  static_cast<unsigned long int>(atoi(argv[1])) );

@#

  /* \newline \aths{  estimate $\EE{R_{i}(n)}$ using |estimate|  \S~\ref{sec:estimateri} } */
  estimate( atof( argv[1]), atof(argv[2]), atof(argv[3]) ) ;

@#

gsl_rng_free( rngtype ); 
return GSL_SUCCESS ; 
}



@* {\bf conclusion and bibliography}.
\label{sec:concl}

The Beta$(\gamma,2-\alpha,\alpha)$-coalescent
\cite{chetwyn-diggle_beta} extends the
Beta$(2-\alpha,\alpha)$-coalescent of \cite{schweinsberg03}; when
$\gamma = 1$ one  recovers  the Beta$(2-\alpha,\alpha$)-coalescent of
\cite{schweinsberg03}.  Moreover, the upper bound affects the
predicted site-frequency spectrum; for any $1 < \alpha < 2$ the
spectrum can be indistinguishable from the one predicted by the
Kingman-coalescent provided $\gamma$ is small enough
\cite{chetwyn-diggle_beta}. The 
Beta$(\gamma,2-\alpha,\alpha)$-coalescent is  determined by
$\gamma \in (0,1]$ and $\alpha \in (1,2)$ and so one should  jointly estimate the
two parameters. The Beta$(\gamma,2-\alpha,\alpha)$-coalescent would be
suitable to compare against population genetic data inherited in a
haploid manner, e.g.\ the site-frequency spectrum of the  mtDNA of diploid populations.





%\bibliographystyle{alpha}
%\bibliography{refs}


\begin{thebibliography}{CDEH25}

\bibitem[CDEH25]{chetwyn-diggle_beta}
JA~Chetwyn-Diggle, Bjarki Eldon
\newblock Beta-coalescents when sample size is large.
\newblock In preparation, 2025+.

\bibitem[DK99]{DK99}
P~Donnelly and T~G Kurtz.
\newblock Particle representations for measure-valued population models.
\newblock {\em Ann Probab}, 27:166--205, 1999.

\bibitem[KL94]{knuth1994cweb}
Donald~Ervin Knuth and Silvio Levy.
\newblock {\em The CWEB system of structured documentation: version 3.0}.
\newblock Addison-Wesley Longman Publishing Co., Inc., Reading, Massachusetts,
  1994.

\bibitem[KR88]{kernighan1988c}
Brian~W Kernighan and Dennis~M Ritchie.
\newblock The {C} programming language, 1988.

\bibitem[Pit99]{P99}
J~Pitman.
\newblock Coalescents with multiple collisions.
\newblock {\em Ann Probab}, 27:1870--1902, 1999.

\bibitem[Sag99]{S99}
S~Sagitov.
\newblock The general coalescent with asynchronous mergers of ancestral lines.
\newblock {\em J Appl Probab}, 36:1116--1125, 1999.

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
