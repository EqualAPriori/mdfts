# Force Field Format
To see a full example, see [example file](#Example_file) below. 

Can also skip to the [new standards](#New_Standards) section to discuss alternative force field formats.

Otherwise, we will break down the sections/elements of a potential, along with [shorthand](#Shorthand)

The main reason for this force field format is to facilitate the definition of *default* parameters, and offer some shorthand to reduce repetitive typing of fields.

The spec is up for modification, and perhaps at some point we should transition to XML, or at least allow for easy conversion from XML to our format and vice versa. I.e. the XML format used by OpenMM allows for successive over-riding of parameter values.

!!! Note "Note"
    the OpenMM XML standard is able to successively define parameters precisely because it uses a convention of *overriding* variables instead of *layering* potentials. I.e it is not allowed (or at least not natural) to layer potentials of the same type on the same pair of atom types in OpenMM (i.e. only one LJ interaction per pair, instead of a superposition of 5)! 

Currently supported potentials:

    bond_harmonic
    pair_ljg
    external_sin
    coulomb_smeared

We use `yaml`, and a force field is essentially a definition with a block for each potential type, and all defaults, parameters, and specific instances of the potential listed said `yaml` block.

## Parts of a potential definition.
Consider the following minimal definition of an external sinusoidal potential, with parts annotated
``` yaml
external_sin:                   # name of the potential-type
    defaults_sim:               # default parameters
        UConst: 0.0 fixed       # following `sim`, this is a fixable (optimizable) parameter
        NPeriods: 1.0 fixed
        PlaneLoc: 0.0           # following `sim`, this is a non-fixable (optimizable) module variable
        PlaneAxis: 0

    params_sim:                 # actual specification
        - SO4 UConst;1.0        # BeadName ParameterName;ParameterValue
```

## Shorthand
Conventions for writing/specifiying parameters.

There are two types of variables in `sim`: optimizable parameters and non-optimizable(fixable) parameters. 

The former are specified by three terms:

    parameter_name, value, fixed/free

While the latter only have two terms:

    parameter_name, value

The second case is a simpler version of the first case, so we only cover the first case.

In fully annotated form, a potential type will consist of a list of dictionaries that fully specify each instance of a potential:

    [
        {   'species': [['btype1','btype2'],['btype1']],
            'name': 'potential_name_typeidentifier',
            'ParamName1': {val:X, fixed:'fixed/free'},
            'ParamName2': {val:Y, fixed:'fixed/free'}
        },
        {
            ...potential 2...
        }
    ]

Or, in python output for a single entry:
``` python
OrderedDict([ ('species', [['SO4'], ['S','C']]), 
                ('name', 'bond_SO4_S'), 
                ('Dist0', {'fixed': True, 'val': 0.0}),
                ('FConst', {'fixed': False, 'val': 1.0})
])
```

The purpose of using dictionaries is that it makes it really transparent to access values of a definition of a potential by `keys`, and it also does not require the user to define every part of a potential in a prescribed order. The dictionary is probably also easiest to make correspondence with `objects` in the future.

The premise is that we want to turn each dictionary entry representing an instance of a potential in a simple format that can even be put on a single line. To do so we introduce several shorthands for making it easier to write the above out in a single line:

- allow for the specification of defaults
- use spaces to distinguish each term of the list
- denote the parameter tuples as a single "word" delimited by ";"'s
- missing values are inferred from the defaults

The user is free to represent things with the shorthand format or any `yaml`-interpretable multi-line form that explicitly denotes lists, dictionaries, etc. In the limit where one uses a fully `yaml`-interpretable specification (with all the `[]{}` brackets), our parser won't be needed, and one can directly work with the resulting `dict-like` object returned by the `yaml` loader.

We start by standardizing an equivalent **`list-form`** representation of each dictionary:

    [ beadtype1, beadtype2, ..., [ParamName,Value,Fixed/Free], [ParamName,Value,FixedFree], ... ]

Note that potential names can be constructed following convention, so we don't have to worry about explicitly including a name/identifier in the above list.

### Representing beadtypes
Let us first address how to denote bead types of a potential, i.e. the dictionary entry key-value pair

    'species': [['SO4'], ['S','C']]

Each potential type knows whether it's a 1-body, 2-body, etc. interaction. An `n`-body potential will correspondingly expect a list of `n` groups/lists/tuples. Since the number of entries/groups is well-defined, we simply specify that the first `n` elements of the **list-form** of the potential contain the bead types. 

Note that we can define groups of more than one bead type, like in the above example, which means to define both a `SO4-S` and a `SO4-C` interaction.

Correspondingly, the above bead type specification for the species involved in this potential can be written in the following equivalent ways that get parsed to the same thing:

    [SO4], [S,C]
    SO4, [S,C]
    [SO4], S;C
    SO4, S;C

where we have used `;` without spaces as a shorthand to concatenate the two elements of the `['S','C']` list. The above entries are meant to to into a '`yaml`' list with `[ .. ]` brackets. Note that `yaml` doesn't need quotes for strings.

Finally, for a completely *inline* representation that is purely a string with no brackets, we can write the above as the beginning of a *inline string* representation:

    SO4 S;C

!!! danger "CONSIDER"
    what do we do if multiple potentials are defined for a pair of beadtypes? Do we override or layer potentials?

    OpenMM, most MD programs, PFTS don't allow for layering potentials, so they are based on overriding potentials. However, sim allows for layering potentials... so it's a different convention.


### Representing parameters
Each parameter is in the most verbose form a key-value pair:

    'Dist0': {'fixed': True, 'val': 0.0}

Where the `value` is itself another dictionary.

We allow the above to be written in several different ways for use within a *`yaml` list*:

    [Dist0, 0.0, fixed]
    {Dist0: 0.0 fixed}
    {Dist0: 0.0 True}
    {Dist0: 0.0;True}

The main shorthand is that when using a dictionary representation for the entry, instead of making the value itself a dictionary with the `fixed` and `val` keys, we allow for `space-` or `,` delimited string representation.

And for a completely inline representation:

    Dist0;0.0;fixed

!!! Note "Note"
    Note that any parameters that aren't defined resort to default values defined in the `defaults_sim` section of each `potential type`.

    Also note that either `fixed/free` or `True/False` can be used.

### Putting it all together
So, the full dictionary

``` python
OrderedDict([ ('species', [['SO4'], ['S','C']]), 
                ('name', 'bond_SO4_S'), 
                ('Dist0', {'fixed': True, 'val': 0.0}),
                ('FConst', {'fixed': False, 'val': 1.0})
])
```

can be written as the following equivalent ways:
    ``` yaml
    [SO4, [S,C], Dist0;0.0;True, FConst;1.0;fixed]          #list-format
    [[SO4], [S,C], [Dist0,0.0,True], {FConst: 1.0 fixed}]   #list-format
    SO4 S;C Dist0;0.0;True FConst;1.0;fixed                 #inline format
    ```

!!! Question "Question"
    When using the `yaml` list to enapsulate the entire entry (so the first line in the above codeblock), am I allowed to use shorthand for the parameters? Or do I have to use lists/dictionaries to specify the parameter?

## New_Standards
STUB

## Example_file
???+ note "Example_file.dat"
    ``` yaml
    # Example yaml force field definition 
    bond_harmonic:
        defaults_sim:
            FConst: 1.0
            Dist0: fixed

        params_sim:
        - SO4 S FConst;1.0;False Dist0,0.0,True
        - C2 S FConst;3.0 Dist0;True
        - S S FConst;100.0 Dist0;1.618
        params_md: #the typical standard

    pair_ljg:
        defaults_sim: #used only by sim
            B: 0.0 True
            Kappa: 1.0 True
            Dist0: 0.0 True
            a_smear: #For now, behavior is to default Kappa, but not pre-given Kappas
                SO4: 0.313173
                S: 0.45
                C2: 0.5

        params_sim: # B*exp(-kappa*(r-Dist0)) + 4 Epsilon ( (Sigma/r)^12 - (Sigma/r)^6 )
            - SO4 S B;1.0;False 
                #- SO4 S;C2 B;2.0 Kappa;0.25
                #- SO4 SO4;S B;2.0
            - SO4 SO4 B;2.0
            - [C2, SO4, [B, 5.0, fixed], {Sigma: 1.0 True}]
            - S S B;1.0
            - C2 C2 B;1.0

    external_sin:
        defaults_sim:
            UConst: 0.0 fixed
            NPeriods: 1.0 fixed
            PlaneLoc: 0.0
            PlaneAxis: 0

        params_sim:
            - SO4 UConst;1.0

    coulomb_smeared: #automatically turns on ewald, then adds the correction potential
        defaults_sim:
            Coef: 1.0 True
            BornA: 1.0 True
            Cut: 1.0
            Shift: True
            ExcludeBondOrd: 0 # only implemented in ewald. currently !=0 values not supported
            a_smear: #first implementation: only overrides if BornA not defined below. Use same smearing convention as Gaussian above, not Born-A convention.
                SO4: 0.5
                S: 0.5
                C2: 0.7

        params_sim: #if BornA = 0, then same as unsmeared Coulomb. Even if no interactions defined, still add ewald.
            - SO4,S SO4,S
            - SO4,S,C2 C2 BornA;2.0

        params_md:
    ```