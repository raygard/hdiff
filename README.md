## hdiff

`hdiff` is a version of diff that implements the _histogram_ algorithm as in jgit.
It duplicates the output of jgit's histogram diff, and usually duplicates (or nearly duplicates) output of `git --histogram`.

It operates in a manner similar to POSIX `diff`, but only supports the default `diff` output and "unified" output formats, and can only compare files (not files to directories or directories to directories), though it does accept `-` for either filename to signify stdin.

`hdiff` does not support other diff options such as `-b` (relaxing whitespace comparison) or `-r` (recursive directory compare), etc.

If enough interest is indicated in the Issues tab, I may consider trying to implement more features.

The code is standard C99, but also needs POSIX extensions to get file modification time in nanoseconds.

There is a fairly complete explanation of the algorithm at the end of the file.

I have also written a [couple](https://raygard.net/2025/01/28/how-histogram-diff-works/) of [posts](https://raygard.net/2025/01/29/a-histogram-diff-implementation/) on this on my [blog](https://raygard.net/).

### Building

```
gcc -O2 -std=c99 -Wall -Wextra -W -Wpointer-arith -Wstrict-prototypes -pedantic -D_POSIX_C_SOURCE=200809L hdiff.c -s -o hdiff
```

Using -O3 works but produces a bogus "warning: ‘kblo’ may be used uninitialized [-Wmaybe-uninitialized]".

### License

Note that the 0BSD license does not require you to keep the license with the code.
Please feel free to copy and re-use the code if you wish.
An attribution to me is nice but not legally required.
