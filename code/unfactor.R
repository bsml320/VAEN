unfactor = function (var) 

{

    if (is.factor(var)) {

        tmp = names(var)

        tmpopt = getOption("warn")

        options(warn = -1)

        out = as.numeric(levels(var))[as.integer(var)]

        options(warn = tmpopt)

        if (any(is.na(out)) & any(is.na(out) != is.na(var))) {

            out = as.character(levels(var))[as.integer(var)]

        }

        names(out) = tmp

    }

    else if (is.data.frame(var)) {

        out = var

        for (i in 1:dim(var)[2]) {

            out[, i] = unfactor(var[, i])

        }

    }

    else if (is.list(var)) {

        out = lapply(var, unfactor)

    }

    else {

        out = var

    }

    return(out)

}
