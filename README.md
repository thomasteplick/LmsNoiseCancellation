<h3>LMS Adaptive Noise Canceller</h3>

<p>
This program is a web application written in Go which uses the html/template package to display the HTML.  To start the server, issue "go run cancellation.go" in the Command Prompt.  To create an adaptive noise canceller using LMS, enter http://127.0.0.1:8080/lmsnoisecanceller in the web browser address bar.    Enter the LMS Equalization Options form data.  The adaptive filter coefficients are averaged over the specified number of trials to produce an ensemble average in order to reduce the variance.  The desired signal is a sine wave in white Gaussian noise.  The user specifies the frequency of the sine wave.  The level of noise is specified by the signal-to-noise ratio (SNR). A reference input is supplied to the adaptive filter input.  The reference input is noise that is correlated with the white noise added to the sine wave.  The correlated noise is produced by applying the white Gaussian noise to a two-pole (complex conjugates) IIR filter, whose pole locations are specified by the user. Poles with radius close to the unit circle produce the most correlation of the white noise and therefore the best filter.
The pole angle does not seem to matter as much as far as removing the noise is concerned.
</p>
<p>
After the adaptive filter is produced, the impulse or frequency response of the adaptive filter can be plotted.  The link <i>Plot Response</i> will take you to the page where you can choose the impulse (time domain) or frequency response.  For the frequency response, the FFT size, number of segments, and window can be chosen.  For this program, one segment with a rectangular window and an FFT size of 8192 is sufficient.  The link <i>Filter Signal</i> will take you to a page where you can display the various stage inputs and outputs. The Filter Signal Options have to be the same as the LMS Equalization Options you used to create the filter.  You have the choice to view the Channel Input, Channel Output, or Noise Canceller Output.  The Channel Input and Channel Output refer to the noisy channel.  Note that each form submittal will have a different input signal sequence. So the sine with random phase will vary from one form submittal to the next. The performance measure for the adaptive filter is how closely the output is to +/- 1 as shown by the Channel Input.
</p>

<h4>LMS Noise Canceller Main Web Page</h4>

![LMSNoiseCanceller1](https://github.com/thomasteplick/LmsNoiseCancellation/assets/117768679/b6ce69f2-7dfd-4ca5-b369-bddd27ad9b74)

<h4>Plot Adaptive Filter Impulse or Frequency Response</h4>

![PlotResponse2](https://github.com/thomasteplick/LmsNoiseCancellation/assets/117768679/09eab446-29b5-4002-ac84-bb335aeb69cd)

<h4>Impulse Response of Adaptive Filter</h4>

![TimeDataPlot3](https://github.com/thomasteplick/LmsNoiseCancellation/assets/117768679/bf48de6b-3d5f-4c2f-a40c-41c7e9778c6b)

<h4>Frequency Respone of Adaptive Filter</h4>

![FrequencyDomainPlot4](https://github.com/thomasteplick/LmsNoiseCancellation/assets/117768679/7177d4bd-bdda-44f3-8392-6464f449c8ba)

<h4>Filter Signal, Channel Input, 10 Hz unit amplitude Sine Wave</h4>

![FilterSignalChannelIn5](https://github.com/thomasteplick/LmsNoiseCancellation/assets/117768679/a446d7f9-7eb1-40cf-b8ae-9c82e2619f51)



