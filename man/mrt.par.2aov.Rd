\name{mrt.par.2aov}
\alias{mrt.par.2aov}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Parametric/non parametric test evaluation
%%  ~~function to do ... ~~
}
\description{
Get info about normality and variances of the given vectors in ordere to decide wether use 2-way ANOVA or Kruskal-Wallis
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mrt.par.2aov(x, y, z, type.of.int = "+" or "*")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a vector of number type (also a dataframe column)
  }
  \item{y}{ a vector of factors
  }
  \item{z}{ a vector of factors
  }
  \item{type.of.int}{"+" to compute the analysis on the two categorial variables only or "*" to compute also the analysis on the interaction
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{Shapiro categorie}{A dataframe containing W and p-value coming from shapiro.test of x in the two different categoric conditions}
  \item{Bartlett}{A dataframe containing k-squared value, df and p-value coming from bartlett.test of x in the two different categoric conditions and, if type.of.int = "*", also from the interaction of the two categorical}
  \item{Shapiro interazione}{optional, only if type.of.int = "*". Return the W and p-value coming from the shapiro test of all the possible interaction of y and z over x}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{ matteo.miotto
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##
mrt.par.2aov(x = rnorm(100), y = sample(rep(c(0,1), 50)), z = sample(rep(c("a", "b", "c", "d"), 25)), type.of.int = "*")

## The function is currently defined as
mrt.par.2aov(x, y, z, type.of.int =)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
