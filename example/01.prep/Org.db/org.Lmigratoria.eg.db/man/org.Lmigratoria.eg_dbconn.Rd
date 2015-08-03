
\name{org.Lmigratoria.eg_dbconn}

\alias{org.Lmigratoria.eg_dbconn}
\alias{org.Lmigratoria.eg_dbfile}
\alias{org.Lmigratoria.eg_dbschema}
\alias{org.Lmigratoria.eg_dbInfo}


\title{Collect information about the package annotation DB}

\description{
  Some convenience functions for getting a connection object to (or collecting
  information about) the package annotation DB.
}

\usage{
  org.Lmigratoria.eg_dbconn()
  org.Lmigratoria.eg_dbfile()
  org.Lmigratoria.eg_dbschema(file="", show.indices=FALSE)
  org.Lmigratoria.eg_dbInfo()
}

\arguments{
  \item{file}{
    A connection, or a character string naming the file to print to (see
    the \code{file} argument of the \code{\link[base]{cat}} function for
    the details).
  }
  \item{show.indices}{
    The CREATE INDEX statements are not shown by default.
    Use \code{show.indices=TRUE} to get them.
  }
}

\details{
  \code{org.Lmigratoria.eg_dbconn} returns a connection object to the
  package annotation DB.  IMPORTANT: Don't call
  \code{\link[DBI:dbConnect]{dbDisconnect}} on the connection object
  returned by \code{org.Lmigratoria.eg_dbconn} or you will break all the
  \code{\link[AnnotationDbi:AnnDbObj-class]{AnnDbObj}} objects defined
  in this package!

  \code{org.Lmigratoria.eg_dbfile} returns the path (character string) to the
  package annotation DB (this is an SQLite file).

  \code{org.Lmigratoria.eg_dbschema} prints the schema definition of the
  package annotation DB.

  \code{org.Lmigratoria.eg_dbInfo} prints other information about the package
  annotation DB.
}

\value{
  \code{org.Lmigratoria.eg_dbconn}: a DBIConnection object representing an
  open connection to the package annotation DB.

  \code{org.Lmigratoria.eg_dbfile}: a character string with the path to the
  package annotation DB.

  \code{org.Lmigratoria.eg_dbschema}: none (invisible \code{NULL}).

  \code{org.Lmigratoria.eg_dbInfo}: none (invisible \code{NULL}).
}

\seealso{
  \code{\link[DBI:dbSendQuery]{dbGetQuery}},
  \code{\link[DBI:dbConnect]{dbConnect}},
  \code{\link[AnnotationDbi:AnnDbObj-class]{dbconn}},
  \code{\link[AnnotationDbi:AnnDbObj-class]{dbfile}},
  \code{\link[AnnotationDbi:AnnDbObj-class]{dbschema}},
  \code{\link[AnnotationDbi:AnnDbObj-class]{dbInfo}}
}

\examples{
  ## Count the number of rows in the "genes" table:
  dbGetQuery(org.Lmigratoria.eg_dbconn(), "SELECT COUNT(*) FROM genes")

  ## The connection object returned by org.Lmigratoria.eg_dbconn() was
  ## created with:
  dbConnect(SQLite(), dbname=org.Lmigratoria.eg_dbfile(), cache_size=64000,
  synchronous=0)

  org.Lmigratoria.eg_dbschema()

  org.Lmigratoria.eg_dbInfo()
}

\keyword{utilities}
\keyword{datasets}

