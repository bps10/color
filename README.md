# COLOR

*These scripts simulate aspects of human color vision based on the Neitz theory of color opponency.*

## Motivation

A parsimonious account of the distribution of monochromatic lights identified as unique, most notably green, has eluded models of color perception.  The standard model in the field consists of three stages predicated upon opponent signals between L- versus M-cones and S-cones versus L- plus M-cone signals.  This model fails to predict the variability in unique hue matching experiments.  We have developed a three stage zone model that takes into account L/M cone ratio to produce distributions of unique hues with accurate statistical characteristics.  The first stage of our model is the known sensitivity functions of the cone photoreceptors.  The second stage is composed of a center of either an L- or M-cone opposed to a surround of all three cones with varying contribution.  In the third stage these opponent signals are weighted based on the probability of occurrence and summated to construct valence curves for the blue-yellow and red-green mechanisms.  Varying the contribution of L- and M-cones to the surround based on an observed distribution of L:M ratios \cite{Carroll2000}, we accurately predict the distributions of unique blue, yellow and green.  Blue and yellow form narrow distributions around 470nm and 578nm, respectively, while green broadly distributes between 495 and 555nm and closely fits the unique green distribution of \cite{Volbrecht1997}. We conclude that our model provides a convincing description of color processing and offers insight into the contribution of L/M ratio in color perception. 

## Color space

1. *Lights* - We adopt the primaries used by WD Wright [r: 650.0nm, g: 530nm, b: 460nm].  As an alternative option we also include those used by Stockman [r: 645nm, g: 526nm, s: 444nm].

2. *Spectral sensitivities* - The default setting uses the Neitz template.  Stockman sensitivities are also available.

3. *Lens & Macula filters* - Stockman's measurements are adopted.

4. *Copunctual points* - LMS = M * rgb, where LMS = [100]' represents protan, LMS = [010]' deutan and LMS = [001]' tritan conditions.

5. *[Protanope, Deuteranope, Tritanope] confusion lines* - first found by solving for the copunctual points.

...

## Color model