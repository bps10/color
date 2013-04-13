# A color model

*These scripts simulate aspects of human color vision based on the Neitz theory of color opponency.*


## Revisiting opponency

####Young-Helmholtz
####Hering
####Hurvich & Jameson
####DeValois & DeValois
####Stockman & Brainard
####Neitz

## Color space

1. *Lights* - We adopt the primaries used by WD Wright [r: 650.0nm, g: 530nm, b: 460nm].  As an alternative option we also include those used by Stockman [r: 645nm, g: 526nm, s: 444nm].

2. *Spectral sensitivities* - The default setting uses the Neitz template.  Stockman sensitivities are also available.

3. *Lens & Macula filters* - Stockman's measurements are adopted.

4. *Copunctual points* - LMS = M * rgb, where LMS = [100]' represents protan, LMS = [010]' deutan and LMS = [001]' tritan conditions.

5. *[Protanope, Deuteranope, Tritanope] confusion lines* - first found by solving for the copunctual points.

...

## Color model