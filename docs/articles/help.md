# Need help?

### Questions

If after repeated tries, checking the [help
files](https://jiscah.github.io/reference/index.md), the
[examples](https://jiscah.github.io/articles/index.md) and
[vignettes](https://jiscah.github.io/articles/vignette_main/book/index.md),
you can’t figure out why sequoia isn’t working the way you think it
should be working: please do feel free to email me at \<jisca.huisman at
gmail.com\>. Usually I’ll get back to you within 1–2 days (during
European daytime).

For example:

- It won’t accept the input data
  - Are you really sure there aren’t any weird symbols in the datafile?
- Assignment rate is (close to) zero, while you are sure there are
  relatives in the dataset
  - Are you really sure the genotyping error rate is as low as you claim
    it to be?
  - Are you really sure
    [`sequoia()`](https://jiscah.github.io/reference/sequoia.md) is able
    to distinguish between different kind of relatives in your dataset?
    Please see
    [here](https://jiscah.github.io/articles/vignette_main/book/faq.md)
    in the main vignette
  - Are you really sure those genotypes belong to relatives? Did you try
    [`CalcPairLL()`](https://jiscah.github.io/reference/CalcPairLL.md) ?
- You are getting an error or warning message you can’t make sense of,
  or which states “Please report bug”
  - Please do email me; something may have slipped through the input
    check. You can check the
    [changelog](https://jiscah.github.io/news/index.md) for bugs that
    have been patched in the development version

### Consultancy work

For several years I worked as freelance consultant, but from September
2025 I have started as researcher at Wageningen University and stopped
with this consultancy work.
