%
%
% BIOINF document class
%
% last revision 2010-02-15 by Ulrich Bodenhofer
%
%
\NeedsTeXFormat{LaTeX2e}[1995/12/01]

\ProvidesClass{bioinf}[2010/02/15]

\def\bioinf@basisclass{article}

\newif\ifbasisclassdefined
\basisclassdefinedfalse

%
% declare options for three basis classes
%
\DeclareOption{article}{\ifbasisclassdefined
\ClassError{bioinf}{You may only supply one basis class.}{Remove the `article' option}
\fi
\def\bioinf@basisclass{article}\basisclassdefinedtrue}

\DeclareOption{report}{\ifbasisclassdefined
\ClassError{bioinf}{You may only supply one basis class.}{Remove the `report' option}
\fi
\def\bioinf@basisclass{report}\basisclassdefinedtrue}

\DeclareOption{book}{\ifbasisclassdefined
\ClassError{bioinf}{You may only supply one basis class.}{Remove the `book' option}
\fi
\def\bioinf@basisclass{book}\basisclassdefinedtrue}

%
% forward all other options to document class article
%
\DeclareOption*{{\ClassError{bioinf}{Option `\CurrentOption' not supported
for BIOINF document class.}{Remove the `\CurrentOption' option}}}

%
% process options
%
\ProcessOptions\relax

%
% load document class article (always with options a4paper and titlepage)
%
\LoadClass[a4paper,11pt,twoside]{\bioinf@basisclass}

\typeout{=============================================================}%
\typeout{BIOINF document class loaded with basis class '\bioinf@basisclass'}%
\typeout{=============================================================}%

%
% font customizations
%
\renewcommand{\encodingdefault}{T1}
\renewcommand{\rmdefault}{ptm}
\renewcommand{\sfdefault}{phv}
\renewcommand{\ttdefault}{cmtt}
\fontencoding{T1}\selectfont

%
% customization of item marks
%
\renewcommand{\labelitemi}{\rule[0.2ex]{0.8ex}{0.8ex}}
\renewcommand{\labelitemii}{$\m@th\bullet$}

%
% include of additional packages
%
\RequirePackage{color}
\RequirePackage{graphicx}
\RequirePackage{ifthen}

%
% setup of page format
%
\setlength{\evensidemargin}{-1in}
\setlength{\oddsidemargin}{-1in}
\setlength{\topmargin}{-1in}
\addtolength\evensidemargin{30mm}
\addtolength\oddsidemargin{30mm}
\addtolength\topmargin{28mm}
\setlength{\textwidth}{150mm}
\setlength{\headheight}{12pt}
\setlength{\headsep}{8mm}
\setlength{\textheight}{222mm}

\ifx\pdfoutput\undefined
\relax
\else
\setlength{\pdfpagewidth}{21cm}
\setlength{\pdfpageheight}{29.7cm}
\fi

%
% setup of customized headings style and listings
%
\input{bioinf-\bioinf@basisclass.txi}

%
% set page style
%
\pagestyle{headings}

%
% adjust parskip
%
\setlength{\parskip}{5pt plus 2pt minus 1pt}

%
% define colors and color commands
%
\definecolor{biowh}{rgb}{1.00,1.00,1.00}
\definecolor{biobl}{rgb}{0.00,0.00,0.00}
\definecolor{bioaz}{rgb}{0.204,0.325,0.631}
\definecolor{biove}{rgb}{0.24,0.63,0.37}
\definecolor{bioli}{rgb}{0.667,0.9333,0.667}

\newcommand{\blue}{\color{bioaz}}
\newcommand{\green}{\color{biove}}

%
% define BIOINF title page commands
%
\newcommand{\customtitlepage}[2]{%
\thispagestyle{empty}
#2
\vfill
\noindent\makebox[0pt][l]{%
\setlength{\unitlength}{1cm}
\begin{picture}(0,0)(2,2.4)
%\put(0,25.5){\makebox(0,0)[lb]{\color{bioli}\rule{19cm}{2.2cm}}}
\put(19,25.5){\makebox(0,0)[rb]{\includegraphics[width=19cm]{bioinf-bar}}}
\put(0.3,26.0){\makebox(0,0)[l]{\color{bioaz}\fontfamily{phv}\bfseries\large%
Institute of Bioinformatics, Johannes Kepler University Linz}}
\put(0.3,26.9){\makebox(0,0)[l]{\color{bioaz}\fontfamily{phv}\bfseries\Huge #1}}
\put(0,1.6){\makebox(0,0)[lb]{\color{bioli}\rule{19cm}{1mm}}}
\put(0,0){\makebox(0,0)[lb]{\fontfamily{phv}\normalsize\begin{tabular}[b]{l}
Institute of Bioinformatics\\
Johannes Kepler University Linz\\
A-4040 Linz, Austria
\end{tabular}}}
\put(19,0){\makebox(0,0)[rb]{\fontfamily{phv}\normalsize\begin{tabular}[b]{r}
Tel. +43 732 2468 4520\\
Fax +43 732 2468 4539\\
http://www.bioinf.jku.at
\end{tabular}}}
%\put(0,24.5){\makebox(0,0)[lt]{\noindent\begin{minipage}{19cm}
%#2
%\end{minipage}}}
\end{picture}}\newpage}

\def\bioinf@number{}
\def\bioinf@keywords{}
\def\bioinf@abstract{}
\def\bioinf@citation{}
\def\bioinf@affiliation{}

\newcommand{\affiliation}[1]{\def\bioinf@affiliation{#1}}
\newcommand{\pubnumber}[1]{\def\bioinf@number{#1}}
\newcommand{\setkeywords}[1]{\def\bioinf@keywords{#1}}
\newcommand{\setabstract}[1]{\def\bioinf@abstract{#1}}
\newcommand{\setcitation}[1]{\def\bioinf@citation{#1}}
\newcommand{\affilmark}[1]{\textsuperscript{\scriptsize #1}}

\newcommand{\doi}[1]{%
\ifx\href\undefined
doi:#1
\else
\href{http://dx.doi.org/#1}{doi:#1}
\fi}

\newcommand{\email}[1]{\ifx\href\undefined
{\em #1}
\else
\href{mailto:#1}{{\em #1}}
\fi}

\newcommand{\reporttitlepage}[1]{\customtitlepage{#1}{%
\mbox{ }
\vfill
\begin{center}
\sffamily\bfseries\LARGE\@title
\end{center}

\begin{center}
{\sffamily\bfseries\large\@author}\\[2mm]
\ifx\@empty\bioinf@affiliation
   \relax
\else
{\sffamily\normalsize\bioinf@affiliation}
\end{center}
\fi

\vspace{1.5cm}

\ifthenelse{\equal{}{\bioinf@abstract}}{%
   \relax
}{%
    \centerline{\rule{15cm}{1pt}}\vspace{2mm}
    \centerline{\begin{minipage}[b]{15cm}
    \noindent{\bf Abstract --- } \bioinf@abstract
    \end{minipage}}
    \ifx\@empty\bioinf@keywords
    \vspace{-1.5mm}\centerline{\rule{15cm}{1pt}}
    \fi
}

\ifthenelse{\equal{}{\bioinf@keywords}}{%
   \relax
}{%
    \ifx\@empty\bioinf@abstract
    \centerline{\rule{15cm}{1pt}}\vspace{2mm}
    \else
    \vspace{5mm}
    \fi
    \centerline{\begin{minipage}[b]{15cm}
    \noindent{\bf Key words --- } \bioinf@keywords
    \end{minipage}}
    \vspace{-1.5mm}\centerline{\rule{15cm}{1pt}}
}

\ifthenelse{\equal{}{\bioinf@citation}}{%
   \relax
}{%
    \vspace{2cm}
    \centerline{\begin{minipage}[b]{15cm}
    \noindent{\bf Citation: } \bioinf@citation
    \end{minipage}}
}

\vfill}}

\newcommand{\researchreporttitlepage}{\reporttitlepage{Research Report \bioinf@number}}

\newcommand{\techreporttitlepage}{\reporttitlepage{Technical Report \bioinf@number}}

\newcommand{\preprinttitlepage}{\reporttitlepage{Publication Preprint}{}}

\newcommand{\reprinttitlepage}{\reporttitlepage{Publication Reprint}{}}

\newcommand{\manualtitlepage}[1][]{\customtitlepage{Software Manual}{%
\mbox{ }
\vfill
\begin{center}
\sffamily\bfseries\LARGE\@title
\end{center}

\ifthenelse{\equal{\@author}{}}{}{\begin{center}
{\sffamily\bfseries\large\@author}\\[2mm]
\ifx\@empty\bioinf@affiliation
   \relax
\else
{\sffamily\normalsize\bioinf@affiliation}
\end{center}
\fi}

\vspace{1cm}

\begin{center}
{\Large\sffamily\bfseries #1}
\end{center}

\vfill

\mbox{ }}}

\newcommand{\lecturenotestitlepage}[1][]{\customtitlepage{Lecture Notes}{%
\ifthenelse{\equal{#1}{}}{%
\vspace*{2cm}
\begin{center}
\sffamily\bfseries\huge \@title
\end{center}

\vspace{2cm}

\begin{center}
\sffamily\bfseries\Large \@author
\end{center}

\vfill

\mbox{ }}{%
{\LARGE\sffamily\bfseries\mbox{ }

\vfill

\begin{center}
#1
\end{center}
\vfill

\vfill

\mbox{ }}}}}

\newcommand{\copyrightpage}[1][]{%
\clearpage
\mbox{ }

\vfill

\vfill

\noindent{\large \textcircled{c}
\ifthenelse{\equal{#1}{}}{%
\@author}{%
#1}}\\[7mm]
\noindent \BIOINFCopyrightNotice\thispagestyle{empty}

\vfill

\mbox{ }
\clearpage}

\newcommand{\BIOINFCopyrightNotice}{This material, no matter whether in printed or electronic
form, may be used for personal and educational use only. Any
reproduction of this manuscript, no matter whether as a whole or in
parts, no matter whether in printed or in electronic form, requires
explicit prior acceptance of the author.}

\endinput
%
% end of scch-tr.cls
%
