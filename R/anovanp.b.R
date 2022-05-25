
anovaNPClass<-if (requireNamespace('jmvcore')) R6::R6Class(
    "anovaNPClass",
    inherit=anovaNPBase,
    private=list(
        #### Member variables ----
        .dataProc=NA,
        .pairsList=NA,
        .kwResults=NA,
        .jtResults=NA,
        .dunnResults=NA,
        .coimResults=NA,
        .dscfResults=NA,
        
        #### Init + run functions ----
        .init=function() {

            private$.dataProc<-NULL
            private$.pairsList<-NULL
            private$.kwResults<-NULL
            private$.jtResults<-NULL
            private$.dunnResults<-NULL
            private$.coimResults<-NULL
            private$.dscfResults<-NULL
            
            deps<-self$options$deps
            group<-self$options$group

            if (is.null(group) || length(deps) == 0)
                return()

            phTables<-self$results$postHoc
            pairs<-self$pairsList

            for (i in seq_along(deps)) {
                pht<-phTables[[i]]

                for (nc in seq_along(pairs))
                    pht$addRow(rowKey=pairs[[nc]],
                               values=list('pair1'=pairs[[nc]][1],
                                           'pair2'=pairs[[nc]][2]))
            }
        },
        .run=function() {

            if(! self$ready)
                return()

            private$.popKwTable()

            private$.popJtTable()
            private$.popDunnPhTables()
            private$.popCoimPhTables()
            private$.popDscfPhTables()

            private$.prepareQQPlot()            

        },
        #### compute functions ----
        .computeKw=function() {

            if ( ! is.null(private$.kwResults) )
                return()

            kwTest<-function(dep, group) {
                if (length(dep) != length(group))
                    jmvcore::reject(
                        jmvcore::format(.("The variables '{}' and '{}' must have the same length."),
                                        dep, group), code='')

                k<-nlevels(group)
                xrank<-rank(dep)
                trank<-table(xrank)
                ni<-tapply(!is.na(dep), group, length)
                N<-sum(ni)
                nr<-length(xrank)
                
                ties<-min(1, 1-sum(trank^3-trank)/(nr^3-nr))
                stat<-(12/(N*(N+1)))*sum(tapply(xrank, group, "sum")^2/ni)-3*(N+1)
                sties<-stat/ties
                pval<-pchisq(stat, df=k-1, lower.tail=FALSE)
                es<-stat*(N+1)/(N^2-1)

                res<-list(stat=as.numeric(stat),
                          sties=as.numeric(sties), 
                          df=as.numeric(k-1), 
                          pval=as.numeric(pval), 
                          es=as.numeric(es))

                return(res)
            }    

            group<-self$options$group
            data<-self$dataProc

            r<-list()
            for (dep in self$options$deps) {
                df<-data.frame(dep=data[[dep]], group=data[[group]])
                if (self$options$miss != 'listwise')
                    df<-na.omit(df)

                r[[dep]]<-kwTest(dep=df$dep, group=df$group)
            }
            return(r)
        },
        .computeJt=function() {

            if ( ! is.null(private$.jtResults) )
                return()

            jtTest<-function(dep, group) {
                if (length(dep) != length(group))
                    jmvcore::reject(
                        jmvcore::format(.("The variables '{}' and '{}' must have the same length."),
                                        dep, group), code='')

                x<-table(as.character(group), dep)
                n_col<-ncol(x)
                n_row<-nrow(x)
                jt_sum<-0

                for(j in 1:(n_col-1))
                    for(i in 1:(n_row-1))
                        jt_sum<-jt_sum+x[i, j]*(0.5*sum(x[(i+1):n_row, j])+sum(x[(i+1):n_row, (j+1):n_col]))

                for(k in 1:(n_row-1))
                    jt_sum<-jt_sum+x[k, n_col]*0.5*sum(x[(k+1):n_row, n_col])

                n<-sum(x)
                nip<-apply(x, 1, sum) 
                npj<-apply(x, 2, sum )
                expect<-(n^2-sum(nip^2))/4

                u1<-n*(n-1)*(2*n+5)-sum(nip*(nip-1)*(2*nip+5))-sum(npj*(npj-1)*(2*npj+5))
                u2<-sum(nip*(nip-1)*(nip-2))*sum(npj*(npj-1)*(npj-2))
                u3<-sum(nip*(nip-1))*sum(npj*(npj-1))

                jt_var<-u1/72+u2/(36*n*(n-1)*(n-2))+u3/(8*n*(n-1))
                jt_se<-sqrt(jt_var)
                jt_z<-(jt_sum-expect)/sqrt(jt_var)

                ## Alternative hypothesis
                pvald<-pnorm(jt_z)                   # decreasing
                pvali<-1-pvald                       # increasing
                pvalt<-2*min(pvald, 1-pvald, 0.5)    # two.sided

                res<-list(stat=as.numeric(jt_sum),
                          seval=as.numeric(jt_se),
                          zval=as.numeric(jt_z),
                          pvalt=as.numeric(pvalt),
                          pvald=as.numeric(pvald),
                          pvali=as.numeric(pvali))

                return(res)
            }    

            group<-self$options$group
            data<-self$dataProc

            r<-list()
            for (dep in self$options$deps) {
                df<-data.frame(dep=data[[dep]], group=data[[group]])
                if (self$options$miss != 'listwise')
                    df<-na.omit(df)

                r[[dep]]<-jtTest(dep=df$dep, group=df$group)
            }
            return(r)
        },        
        .computeDunn=function() {

            if ( ! is.null(private$.dunnResults) )
                return()

            dunnPhTest<-function(dep, group) {
                if (length(dep) != length(group))
                    jmvcore::reject(
                        jmvcore::format(.("The variables '{}' and '{}' must have the same length."),
                                        dep, group), code='')

                k<-nlevels(group)
                xrank<-rank(dep)
                trank<-table(xrank)
                ni<-tapply(!is.na(dep), group, length)
                N<-sum(ni)
                rb<-tapply(xrank, group, mean, na.rm=T)
                ties<-sum(trank^3-trank)/(12*(length(xrank)-1))

                ## function for pairwise comparisons
                paircomp<-function(i, j) {
                    dif<-abs(rb[i]-rb[j])
                    A<-N*(N+1)/12
                    B<-(1/ni[i]+1/ni[j])
                    zval<-dif/sqrt((A-ties)*B)
                    return(zval)
                }
                stat<-pairwise.table(paircomp, levels(group), p.adjust.method="none")
                pval<-pnorm(abs(stat), lower.tail=FALSE)
                stat=t(stat)[upper.tri(t(stat), diag=TRUE)]
                pval=t(pval)[upper.tri(t(pval), diag=TRUE)]

                ## p adjust method
                ncomp<-k*(k-1)/2
                pnone<-pval
                pbonf<-pmin(1, pval*ncomp)
                psidak<-pmin(1, 1-(1-pval)^ncomp)
                ## holm method 
                Psort<-matrix(c(pval, 1:ncomp, rep(0, ncomp)), 3, ncomp, byrow=TRUE)
                Psort<-Psort[,order(Psort[1,])]
                for (i in 1:ncomp) {
                    adjust<-ncomp+1-i
                    Psort[1,i]<-pmin(1, Psort[1,i]*adjust)
                    Psort[3,i]<-Psort[1,i] <= 0.025 # alpha/2
                }
                Psort<-Psort[,order(Psort[2,])]
                pholm<-Psort[1,]

                res<-list(stat=as.numeric(stat),
                          pval=as.numeric(pval),
                          pnone=as.numeric(pnone), 
                          pbonf=as.numeric(pbonf),
                          psidak=as.numeric(psidak),
                          pholm=as.numeric(pholm))

                return(res)
            }

            group<-self$options$group
            data<-self$dataProc

            r<-list()
            for (dep in self$options$deps) {
                df<-data.frame(dep=data[[dep]], group=data[[group]])
                if (self$options$miss != 'listwise')
                    df<-na.omit(df)

                r[[dep]]<-dunnPhTest(dep=df$dep, group=df$group)
            }

            return(r)
        },
        .computeCoim=function() {

            if ( ! is.null(private$.coimResults) )
                return()

            coimPhTest<-function(dep, group) {
                if (length(dep) != length(group))
                    jmvcore::reject(
                        jmvcore::format(.("The variables '{}' and '{}' must have the same length."),
                                        dep, group), code='')

                k<-nlevels(group)
                xrank<-rank(dep)
                trank<-table(xrank)
                ni<-tapply(!is.na(dep), group, length)
                N<-sum(ni)
                rb<-tapply(xrank, group, mean, na.rm=T)
                nr<-length(xrank)

                ties<-min(1, 1-sum(trank^3-trank)/(nr^3-nr))
                H<-(12/(N*(N+1)))*sum(tapply(xrank, group, "sum")^2/ni)-3*(N+1)
                adj<-H/ties                

                if (ties == 1)
                    S2<-N*(N+1)/12
                else 
                    S2<-(1/(N-1))*(sum(xrank^2)-(N*(((N+1)^2)/4)))

                ## function for pairwise comparisons
                paircomp<-function(i, j) {
                    dif<-rb[i]-rb[j]
                    B<-(1/ni[i]+1/ni[j])
                    D<-(N-1-adj)/(N-k)
                    tval<-dif/sqrt(S2*B*D)
                    return(tval)
                }
                stat<-pairwise.table(paircomp, levels(group), p.adjust.method="none")
                pval<-pt(q=abs(stat), df=N-k, lower.tail=FALSE)
                stat=t(stat)[upper.tri(t(stat), diag=TRUE)]
                pval=t(pval)[upper.tri(t(pval), diag=TRUE)]

                ## p adjust method 
                ncomp<-k*(k-1)/2
                pnone<-pval
                pbonf<-pmin(1, pval*ncomp)
                psidak<-pmin(1, 1-(1-pval)^ncomp)
                ## holm method 
                Psort<-matrix(c(pval, 1:ncomp, rep(0, ncomp)), 3, ncomp, byrow=TRUE)
                Psort<-Psort[,order(Psort[1,])]
                for (i in 1:ncomp) {
                    adjust<-ncomp+1-i
                    Psort[1,i]<-pmin(1, Psort[1,i]*adjust)
                    Psort[3,i]<-Psort[1,i] <= 0.025 # alpha/2
                }
                Psort<-Psort[,order(Psort[2,])]
                pholm<-Psort[1,]

                res<-list(stat=as.numeric(stat),
                          pval=as.numeric(pval),
                          pnone=as.numeric(pnone), 
                          pbonf=as.numeric(pbonf),
                          psidak=as.numeric(psidak),
                          pholm=as.numeric(pholm))

                return(res)
            }

            group<-self$options$group
            data<-self$dataProc

            r<-list()
            for (dep in self$options$deps) {
                df<-data.frame(dep=data[[dep]], group=data[[group]])
                if (self$options$miss != 'listwise')
                    df<-na.omit(df)

                r[[dep]]<-coimPhTest(dep=df$dep, group=df$group)
            }

            return(r)
        },
        .computeDscf=function() {

            if ( ! is.null(private$.dscfResults) )
                return()

            dscfPhTest<-function(dep, group) {
                if (length(dep) != length(group))
                    jmvcore::reject(
                        jmvcore::format(.("The variables '{}' and '{}' must have the same length."),
                                        dep, group), code='')

                k<-nlevels(group)
                N<-tapply(dep, group, length)
                glev<-levels(group)

                ## function for pairwise comparisons
                paircomp<-function(i, j) {
                    nn<-N[i]
                    m<-N[j]
                    xraw<-c(dep[group==glev[i]], dep[group==glev[j]])
                    rankx<-rank(xraw)
                    lev<-c(group[group==glev[i]], group[group==glev[j]])
                    id<-!glev %in% c(glev[i], glev[j])
                    lev<-droplevels(lev, exclude=glev[id])
                    R<-tapply(rankx, lev, sum)
                    U<-c(m*nn+(m*(m+1)/2), m*nn+(nn*(nn+1)/2))-R
                    t<-table(rankx)
                    V<-(m*nn/((m+nn)*(m+nn-1)))*(((m+nn)^3-(m+nn))/12-sum((t^3-t)/12))
                    q<-sqrt(2)*(min(U)-m*nn/2)/sqrt(V)
                    return(q)
                }

                stat<-pairwise.table(paircomp, levels(group), p.adjust.method="none")
                pval<-ptukey(abs(stat), nmeans=k, df=Inf, lower.tail=FALSE)
                stat=t(stat)[upper.tri(t(stat), diag=TRUE)]
                pval=t(pval)[upper.tri(t(pval), diag=TRUE)]
                res<-list(stat=as.numeric(stat), pval=as.numeric(pval))

                return(res)
            }

            group<-self$options$group
            data<-self$dataProc

            r<-list()
            for (dep in self$options$deps) {
                df<-data.frame(dep=data[[dep]], group=data[[group]])
                if (self$options$miss != 'listwise')
                    df<-na.omit(df)

                r[[dep]]<-dscfPhTest(dep=df$dep, group=df$group)
            }

            return(r)
        },
        #### pop tables functions ----
        .popKwTable=function() {
            deps<-self$options$deps
            kwtable<-self$results$kwTable
            kwRes<-self$kwResults

            rows<-list()
            for (dep in self$options$deps) {
                rows[[dep]]<-list(
                    stat=kwRes[[dep]]$stat,
                    sties=kwRes[[dep]]$sties,
                    df=kwRes[[dep]]$df,
                    pval=kwRes[[dep]]$pval,
                    es=kwRes[[dep]]$es)

                kwtable$setRow(rowKey=dep, rows[[dep]])
            }
        },
        .popJtTable=function() {

            if (! self$options$jttest)
                return()

            deps<-self$options$deps
            jttable<-self$results$jtTable
            jtRes<-self$jtResults

            rows<-list()
            for (dep in self$options$deps) {
                rows[[dep]]<-list(
                    stat=jtRes[[dep]]$stat,
                    seval=jtRes[[dep]]$seval,
                    zval=jtRes[[dep]]$zval,
                    pval=jtRes[[dep]]$pvalt,
                    pvald=jtRes[[dep]]$pvald,
                    pvali=jtRes[[dep]]$pvali)

                jttable$setRow(rowKey=dep, rows[[dep]])
            }
        },
        .popDunnPhTables=function() {

            if (! self$options$dunnpairs)
                return()

            deps<-self$options$deps
            group<-self$options$group

            phTables<-self$results$postHoc
            dunnRes<-self$dunnResults

            pairs<-self$pairsList

            for (i in seq_along(deps)) {

                pht<-phTables[[i]]
                r<-dunnRes[[deps[i]]]

                for (nc in seq_along(pairs)) {
                    if (pht$getCell(rowKey=pairs[[nc]], 'stat[dunn]')$isEmpty) {
                        pht$setStatus('running')
                        private$.checkpoint()

                        pht$setRow(rowKey=pairs[[nc]],
                                   list('stat[dunn]'=r$stat[[nc]],
                                        'pval[dunn]'=r$pval[[nc]],
                                        'pnone[dunn]'=r$pnone[[nc]],
                                        'pbonf[dunn]'=r$pbonf[[nc]],
                                        'psidak[dunn]'=r$psidak[[nc]],
                                        'pholm[dunn]'=r$pholm[[nc]]))

                        pht$setStatus('complete')
                    }
                }
            }
        },
        .popCoimPhTables=function() {

            if (! self$options$coimpairs)
                return()

            deps<-self$options$deps
            group<-self$options$group

            phTables<-self$results$postHoc
            coimRes<-self$coimResults

            pairs<-self$pairsList

            for (i in seq_along(deps)) {

                pht<-phTables[[i]]
                r<-coimRes[[deps[i]]]

                for (nc in seq_along(pairs)) {
                    if (pht$getCell(rowKey=pairs[[nc]], 'stat[coim]')$isEmpty) {
                        pht$setStatus('running')
                        private$.checkpoint()

                        pht$setRow(rowKey=pairs[[nc]],
                                   list('stat[coim]'=r$stat[[nc]],
                                        'pval[coim]'=r$pval[[nc]],
                                        'pnone[coim]'=r$pnone[[nc]],
                                        'pbonf[coim]'=r$pbonf[[nc]],
                                        'psidak[coim]'=r$psidak[[nc]],
                                        'pholm[coim]'=r$pholm[[nc]]))

                        pht$setStatus('complete')
                    }
                }
            }
        },
        .popDscfPhTables=function() {

            if (! self$options$dscfpairs)
                return()

            deps<-self$options$deps
            group<-self$options$group

            phTables<-self$results$postHoc
            dscfRes<-self$dscfResults

            pairs<-self$pairsList

            fnA <-.('P value adjustment method: single-step')
            fnB <-.('This correction is unnecessary for this test.')

            for (i in seq_along(deps)) {

                pht<-phTables[[i]]
                r<-dscfRes[[deps[i]]]

                for (nc in seq_along(pairs)) {
                    if (pht$getCell(rowKey=pairs[[nc]], 'stat[dscf]')$isEmpty) {
                        pht$setStatus('running')
                        private$.checkpoint()

                        pht$setRow(rowKey=pairs[[nc]],
                                   list('stat[dscf]'=r$stat[[nc]],
                                        'pval[dscf]'=r$pval[[nc]],
                                        'pnone[dscf]'=NaN,
                                        'pbonf[dscf]'=NaN,
                                        'psidak[dscf]'=NaN,
                                        'pholm[dscf]'=NaN))

                        pht$addFootnote(rowKey=pairs[[nc]], col='pval[dscf]', fnA)
                        pht$addFootnote(rowKey=pairs[[nc]], col='pnone[dscf]', fnB)
                        pht$addFootnote(rowKey=pairs[[nc]], col='pbonf[dscf]', fnB)
                        pht$addFootnote(rowKey=pairs[[nc]], col='psidak[dscf]', fnB)
                        pht$addFootnote(rowKey=pairs[[nc]], col='pholm[dscf]', fnB)

                        pht$setStatus('complete')
                    }
                }
            }
        },
        .prepareQQPlot=function() {
            data<-self$dataProc
            plots<-self$results$plots
            group<-self$options$group

            for (dep in self$options$deps) {
                image<-plots$get(key=dep)$qq
                df<-data.frame(dep=data[[dep]], group=data[[group]])
                image$setState(df)
            }
        },
        .qq=function(image, ggtheme, theme, ...) {

            if (is.null(image$state))
                return(FALSE)
 
            # Perform QQ plots by group
            p<-ggplot2::ggplot(data=image$state, 
                               mapping=ggplot2::aes(sample=dep, color=group)) +
                ggplot2::geom_qq() +
                ggplot2::geom_qq_line() +
                ggplot2::facet_wrap(~ group, scales="free") + ggtheme +
                ggplot2::labs(x=.("Theoretical Quantiles"), 
                              y=.("Sample Quantiles")) +
                ggplot2::theme(legend.position="None")

            return(p)
        },
        .cleanData=function() {
            deps<-self$options$deps
            group<-self$options$group
            varNames<-c(deps, group)

            data<-jmvcore::select(self$data, varNames)

            for (dep in deps)
                data[[dep]]<-jmvcore::toNumeric(data[[dep]])

            data[[group]]<-droplevels(as.factor(data[[group]]))

            if (any(deps == group))
                jmvcore::reject(
                    jmvcore::format(.("The grouping variable '{}' must not also be a dependent variable."),
                                    group), code='')

            # exclude rows with missings in the grouping variable
            data<-data[ ! is.na(data[[group]]),]

            lvls<-base::levels(data[[group]])
            if (length(lvls) == 0) {
                jmvcore::reject(
                    jmvcore::format(.("The grouping variable '{}' contains no data."),
                                    group), code='')
            } else if (length(lvls) == 1) {
                jmvcore::reject(
                    jmvcore::format(.("The grouping variable '{}' all observations are in the same group."),
                                    group), code='')
            } else if (length(lvls) == 2) {
                jmvcore::reject(
                    jmvcore::format(.("The grouping variable '{}' should have at least 3 or more levels."),
                                    group), code='')
            }

            if (self$options$miss == "listwise") {
                data<-naOmit(data)
                if (dim(data)[1] == 0)
                    jmvcore::reject(
                        jmvcore::format(.("Grouping variable '{}' has less than 3 levels after missing values are excluded"),
                                        group), code='')
            }
            return(data)
        },
        .genPairs=function() {
            group<-self$options$group
            data<-self$data

            lvls<-base::levels(data[[group]])
            pairsList<-list()
            if (length(lvls) > 0) {
                pairsMatrix<-utils::combn(lvls, 2)
                for (i in seq_len(dim(pairsMatrix)[2]))
                    pairsList[[i]]<-pairsMatrix[,i]
            }
            return(pairsList)
        },
        .sourcifyOption=function(option) {
            if (option$name %in% c('deps', 'group'))
                return('')
            super$.sourcifyOption(option)
        },
        .formula=function() {
            jmvcore:::composeFormula(self$options$deps, self$options$group)
        }
    ), 
    #### Active bindings ----
    active=list(
        ready=function() {
            group<-self$options$group
            deps<-self$options$deps

            return( ! is.null(group)
                    && length(deps) > 0
                    && nrow(self$data) > 0)
        },
        dataProc=function() {
            if (is.null(private$.dataProc))
                private$.dataProc<-private$.cleanData()
            
            return(private$.dataProc)
        },
        pairsList=function() {
            if (is.null(private$.pairsList))
                private$.pairsList<-private$.genPairs()
            
            return(private$.pairsList)
        },
        kwResults=function() {
            if (is.null(private$.kwResults))
                private$.kwResults<-private$.computeKw()
            
            return(private$.kwResults)
        },
        jtResults=function() {
            if (is.null(private$.jtResults))
                private$.jtResults<-private$.computeJt()
            
            return(private$.jtResults)
        },
        dunnResults=function() {
            if (is.null(private$.dunnResults))
                private$.dunnResults<-private$.computeDunn()
            
            return(private$.dunnResults)
        },
        coimResults=function() {
            if (is.null(private$.coimResults))
                private$.coimResults<-private$.computeCoim()
            
            return(private$.coimResults)
        },
        dscfResults=function() {
            if (is.null(private$.dscfResults))
                private$.dscfResults<-private$.computeDscf()
            
            return(private$.dscfResults)
        }
    )
)
