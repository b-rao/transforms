# transforms

A program that can convolve two wav files, the input and impulse together to get an output that sounds like the input being played in the impulse environment.

This was a assignment on optimization. The first try would have you use a very slow time domain convolution,

but then once you implemented a faster algorithm, you can see the improvements.

## about

Tested on linux

Takes two files, the input wav and the impulse response wav, and a specified output filename.

```
./convolve input.wav ir.wav out.wav
```

stereo wav files are not supported.
