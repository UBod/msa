%
%
% UBMANUAL document class
%
% last revision 2024-04-19 by Ulrich Bodenhofer
%
%
\NeedsTeXFormat{LaTeX2e}[1995/12/01]

\ProvidesClass{ubmanual}[2024/04/19]

%
% load document class article (always with options a4paper and titlepage)
%
\LoadClass[a4paper,11pt,twoside]{article}

\typeout{=============================================================}%
\typeout{UBMANUAL document class loaded}%
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
\def\ps@headings{%
      \let\@oddfoot\@empty\let\@evenfoot\@empty
      \def\@evenhead{\underline{\makebox[\textwidth][l]{%
          \sffamily\thepage\hfill\leftmark}}}%
      \def\@oddhead{\underline{\makebox[\textwidth][l]{%
          \sffamily\rightmark\hfill\thepage}}}%
      \let\@mkboth\markboth
    \def\sectionmark##1{%
      \markboth {
        \ifnum \c@secnumdepth >\z@
          \thesection\quad
        \fi
        ##1}{
        \ifnum \c@secnumdepth >\z@
          \thesection\quad
        \fi
        ##1}}%
        \def\subsectionmark##1{}}

\renewenvironment{abstract}{\section*{\abstractname}%
    \@mkboth{\abstractname}{\abstractname}}{}
\newenvironment{acknowledgments}{\section*{Acknowledgments}%
    \@mkboth{Acknowledgments}{Acknowledgments}}{}
\renewcommand\tableofcontents{%
    \section*{\contentsname
        \@mkboth{\contentsname}{\contentsname}}%
    \@starttoc{toc}%
    }
\renewcommand\listoffigures{%
    \section*{\listfigurename}%
      \@mkboth{\listfigurename}%
              {\listfigurename}%
    \@starttoc{lof}%
    }
\renewcommand\listoftables{%
    \section*{\listtablename}%
      \@mkboth{\listtablename}%
         {\listtablename}%
    \@starttoc{lot}%
    }
\renewenvironment{thebibliography}[1]
     {\section*{\refname}%
      \@mkboth{\refname}{\refname}%
      \list{\@biblabel{\@arabic\c@enumiv}}%
           {\settowidth\labelwidth{\@biblabel{#1}}%
            \leftmargin\labelwidth
            \advance\leftmargin\labelsep
            \@openbib@code
            \usecounter{enumiv}%
            \let\p@enumiv\@empty
            \renewcommand\theenumiv{\@arabic\c@enumiv}}%
      \sloppy
      \clubpenalty4000
      \@clubpenalty \clubpenalty
      \widowpenalty4000%
      \sfcode`\.\@m}
     {\def\@noitemerr
       {\@latex@warning{Empty `thebibliography' environment}}%
      \endlist}

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
\definecolor{ubwh}{rgb}{1.00,1.00,1.00}
\definecolor{ubbl}{rgb}{0.00,0.00,0.00}
\definecolor{ubdg}{rgb}{0.30,0.30,0.30}
\definecolor{ubgr}{rgb}{0.55,0.55,0.55}
\definecolor{ublg}{rgb}{0.80,0.80,0.80}

\newcommand{\lightgrey}{\color{ublg}}
\newcommand{\darkgrey}{\color{ubdg}}
\newcommand{\grey}{\color{ubgr}}

%
% define UBMANUAL title page commands
%
\newcommand{\customtitlepage}[3]{%
\thispagestyle{empty}
#2
\vfill
\noindent\makebox[0pt][l]{%
\setlength{\unitlength}{1cm}
\begin{picture}(0,0)(2,2.4)
\put(0,25.5){\makebox(0,0)[lb]{\color{ublg}\rule{19cm}{2.2cm}}}
\put(0.3,26.5){\makebox(0,0)[l]{\fontfamily{phv}\bfseries\Huge #1}}
\put(0,1.6){\makebox(0,0)[lb]{\color{ublg}\rule{19cm}{1mm}}}
\put(0,0){\makebox(0,0)[lb]{\fontfamily{phv}\normalsize\begin{tabular}[b]{l}
#3
\end{tabular}}}
\end{picture}}\newpage}

\def\ubmanual@number{}
\def\ubmanual@keywords{}
\def\ubmanual@abstract{}
\def\ubmanual@citation{}
\def\ubmanual@affiliation{}

\newcommand{\affiliation}[1]{\def\ubmanual@affiliation{#1}}
\newcommand{\affilmark}[1]{\textsuperscript{\scriptsize #1}}

\newcommand{\doi}[1]{%
\ifx\href\undefined
doi:#1
\else
\href{http://doi.org/#1}{doi:#1}
\fi}

\newcommand{\email}[1]{\ifx\href\undefined
{\em #1}
\else
\href{mailto:#1}{{\em #1}}
\fi}

\newcommand{\manualtitlepage}[2]{\customtitlepage{User Manual / Package Vignette}{%
\mbox{ }
\vfill
\begin{center}
\sffamily\bfseries\LARGE\@title
\end{center}

\ifthenelse{\equal{\@author}{}}{}{\begin{center}
{\sffamily\bfseries\large\@author}\\[2mm]
\ifx\@empty\ubmanual@affiliation
   \relax
\else
{\sffamily\normalsize\ubmanual@affiliation}
\end{center}
\fi}

\vspace{1cm}

\begin{center}
{\Large\sffamily\bfseries #1}
\end{center}

\vfill

\mbox{ }}{%
Link to package on GitHub: \href{#2}{#2} \\
\mbox{ } \\
\mbox{ }
}}

\newcommand{\copyrightpage}[1][]{%
\clearpage
\mbox{ }

\vfill

\vfill

\noindent{\large \textcircled{c}
\ifthenelse{\equal{#1}{}}{%
\@author}{%
#1}}\\[7mm]
\noindent \UBMANUALCopyrightNotice\thispagestyle{empty}

\vfill

\mbox{ }
\clearpage}

\newcommand{\UBMANUALCopyrightNotice}{This material, no matter whether in printed or electronic
form, may be used for personal and educational use only. Any
reproduction of this manuscript, no matter whether as a whole or in
parts, no matter whether in printed or in electronic form, requires
explicit prior acceptance of the author.}

\endinput
%
% end of ubmanual.cls
%
