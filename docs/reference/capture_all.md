# Silently evaluate an expression while returning its value

Evaluates an expression, suppressing all console output (print, cat,
implicit printing) and all messages, while still returning the normal
value of the expression.

## Usage

``` r
capture_all(exprs, file = NULL)
```

## Arguments

- exprs:

  Expression to be evaluated.

- file:

  Optional path to a file where captured output/messages should be
  written. If `NULL`, they are discarded to the OS null device.

## Value

The normal return value of `exprs`.
