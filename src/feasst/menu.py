"""
Todo:

Undo seems to have some kind of memory leak

Add src/feasst/feasst.cpp to pydoc

Remember: site type names, site names

Configuration [param=epsilon,sigma,etc] arguments need to be input, with list of defaults plus custom ?

Particle builder: remember names of particles, sites, etc. Obtain sites/bonds/etc for future reference. txt->json particle file.

Print the version with MonteCarlo so exectuable can compare against it (not at constructor, what about prefetch?).
Print MonteCarlo at first line for reproducibility, like text interface script

Or... make more use of the actual program errors to run the simulation in client-server model interactively (undo doesn't work so well in that case)
|- Advise user to add a chemical potential of type if TrialAdd or TrialGrow
|Make user provide a chemical potential (if not xyz or add_num?)

Enable Let, if, for and somehow remember the variables for user input?

For variables, enable them as an option (like defaults but before auto chooser) for all classes. The description is the entire let command.

For loops define temporary variables that go away at endFor (allow nesting of stored for)

For and if can have automatic 2-spacs indenting

See this example for ifs where question marks are needed
https://pages.nist.gov/feasst/plugin/gibbs/tutorial/tutorial_2_trappe_gibbs.html

Document (do all this together with each header):
- parameters c++ comments should have default or optional (as implemented in analyze_public_interface) to distinguish
For HTML doc
- Link to and or repeat information in arguments of one class that are described in another class (e.g., group).
- parse other links
- derived class arguments
- captialize the first character on argument descriptions
- public interface classes need to be one per file to simplify header file analysis
- footcite, biblio... (what if can automate refs and latex and put that in waiter?)
- instructions for documentation. optional/default, This class is [experimental/deprecated], etc. optional arg in base class is overriden by required arg in derived (TrialAdd::particle_type)
"""

import copy
import os
import argparse
import json
import curses
from curses.textpad import Textbox, rectangle
import logging as log
import textwrap
from importlib.resources import files
import subprocess
import webbrowser

# Configure log to a file
#log.basicConfig(filename="feasst-menu.debug.txt", level=log.DEBUG, format='[%(filename)s:%(lineno)d] %(message)s')

# Custom options
LOAD_SAVE_FILE = '[Load save]'
LOAD_SCRIPT_FILE = '[Load script]'
LOAD_PARTICLE_FILE = '[Load particle file]'
OPEN_HTML = '[Open HTML]'
SAVE_FILE = '[Save]'
NO_MORE_ARGUMENTS = '[Finish \"'
NO_MORE_PARTICLES = '[No more particles]'
NUMBER_ARG_ROWS= "[Number of rows for arguments]"
NUMBER_SCRIPT_ROWS= "[Number of rows for script]"
TOGGLE_TUTORIAL= "[Toggle tutorial mode]"
UNIQUE_START = '__start__'
UNIQUE_PARTICLE = '__particle__'
OPTIONAL = '__optional__'
COMMENT_LINE = '[Add comment]'
CUSTOM_LINE = '[Add line]'
RUN_COMMAND = '[Save and run]'

def populate_base(data):
    for cl in data:
        for idx,lst in enumerate(data[cl]):
            if len(lst) > 0:
                cl2 = lst[0]
                #if cl == 'Potential'+OPTIONAL:
                #    print('cl2', cl2, 'dat---', data[cl])
                if cl2 in data and lst[2] == "":  # description is empty so its a class (analyze line ~47)
                    data[cl].pop(idx)
                    if cl[-12:] == OPTIONAL:
                        data[cl] += data[cl2+OPTIONAL]
                    else:
                        data[cl] += data[cl2]
    return data

def gen_data(descript, particles, factory, tutorial):
    # remove deprecated and experimental classes from both descript and factory
    #for cl in reversed(descript):
    #print('*******')
    excludes = ['This class is deprecated', 'This class is experimental', 'This class is for developers', 'This class is in development']
    for icl,cl in reversed(list(enumerate(descript))):
        if cl != "":
            if any(excl in descript[cl]['descript'] for excl in excludes):
            #if 'This class is deprecated' in descript[cl]['descript'] or \
            #   'This class is experimental' in descript[cl]['descript']:
                descript.pop(cl)
                #print(cl)
                for fac in factory:
                    for ider,der in reversed(list(enumerate(factory[fac]))):
                        if der == cl:
                            factory[fac].pop(ider)
                            #print('poping from fac', ider, der, factory[fac])
    #print('*******')
    #print(factory['Action'][13])
    #quit()

    # required and optional arguments
    data = dict()
    for cl in descript:
        if cl != "":
            data[cl] = []
            for _, arg in enumerate(descript[cl]['required']):
                if arg[0] in factory and arg[1] != "":
                    data[cl].append([arg[0], cl, arg[1], ''])
                else:
                    data[cl].append([arg[0], cl, arg[1], ''])
            clo = cl+OPTIONAL
            data[clo] = []
            for _, arg in enumerate(descript[cl]['optional']):
                if arg[1] == "":
                    data[clo].append([arg[0], clo, arg[1], ''])
                    if arg[0] in descript:
                    #if arg[1] == "": # arg[0] means its a class, not argument. Arg[1] means no description (analyze_public_interface.py line ~47), so its not asking for base class as argument.
                        # all base classes end up in optional, so add to required as well
                        data[cl].append([arg[0], cl, arg[1], ''])
                else:
                    data[clo].append([arg[0], clo, arg[1], ''])

    # add arguments from base classes. May need to repeat multiple times?
    for _ in range(5):
        log.debug(str(data['Metropolis']))
        data = populate_base(data)
    log.debug(str(data['Metropolis']))

    # remove optional argument if it also shows up as required (e.g., TrialAdd requires particle_type, which is optional in TrialSelect)
    for cl in data:
        if OPTIONAL not in cl:
            for idx,lst in enumerate(data[cl]):
                for idx2,lst2 in reversed(list(enumerate(data[cl+OPTIONAL]))):
                    if lst[0] == lst2[0]:
                        data[cl+OPTIONAL].pop(idx2)

    # Add default "no more arguments" option (after base population or there will be more than one)
    for cl in descript:
        log.debug('cln:'+cl)
        if cl != "":
            data[cl+OPTIONAL].insert(0, (NO_MORE_ARGUMENTS+cl+"\"]", cl+OPTIONAL, ""))

    # customize particle file selection when Configuration is selected
    class_config = (["Configuration", "particles_select0", descript['Configuration']['descript']])
    cust_particle = (LOAD_PARTICLE_FILE, 'particles_select1', 'Enter the full path and full name of a customized particle file.', '')
    data['particles_select0'] = []
#    cust_particle]
    data['particles_select1'] = [(NO_MORE_PARTICLES, 'Configuration', 'Do not add any more particles to the Configuration, and continue with other optional parameters.')]
#    , cust_particle]
    for part in particles:
        firstline = '# This file is located in /feasst/'+particles[part][1]+'\n'
        data['particles_select0'].append((part, 'particles_select1', firstline+particles[part][0], particles[part][1]))
        data['particles_select1'].append((part, 'particles_select1', firstline+particles[part][0], particles[part][1]))
    data['particles_select0'].append(cust_particle)
    data['particles_select1'].append(cust_particle)
    # remove particle_type
    for idat,dat in enumerate(data['Configuration']):
        if dat[0] == 'particle_type':
            data['Configuration'].pop(idat)
            break;
    # point non-optional config parameters to the optional ones
    for idat,dat in enumerate(data['Configuration']):
        data['Configuration'][idat][3] = 'Configuration' + OPTIONAL

    # build factory options
    for fac in factory:
        data[fac] = []
        for deriv in factory[fac]:
            if deriv in descript:
                hi=1
                data[fac].append((deriv, deriv, descript[deriv]['descript']))

    # search arguments for factory and replace links
    for cl in data:
        for idx,lst in enumerate(data[cl]):
            for fac in factory:#fidx,fac in enumerate(factory):
                if lst[0] == fac:
                    data[cl][idx][1] = fac
                    data[cl][idx][3] = 'factory'
                    #print('lst', lst[0], 'cl', cl, data[cl][idx])

    # redirect arguments to select a list of available options
    for cl in data:
        for idx,lst in enumerate(data[cl]):
            if len(lst) > 0:
                if lst[0] == 'particle_type':
                    data[cl][idx][3] = 'particle_type_chooser'
                if lst[0] == 'config':
                    data[cl][idx][3] = 'configuration_chooser'
                if lst[0] == 'group':
                    if cl == 'Configuration'+OPTIONAL:
                        data[cl][idx][1] = 'Group'
                    else:
                        data[cl][idx][3] = 'group_chooser'

    # expand classes
    class_ = dict()
    for cl in descript:
        class_[cl] = ([cl, cl, descript[cl]['descript']])
    for fac in factory:
        class_[fac] = ([fac, fac, descript[fac]['descript']])

    every_mc_no_for = [
      class_['Checkpoint'],
      ([COMMENT_LINE, "", "Add a comment line to improve the readability of the script. Comment lines that begin with a \"#\" symbol are ignored by the FEASST executable.", ""]),
      ([CUSTOM_LINE, "", "Enter a line directly into the script that will be read by the FEASST executable (not recommended, although this can also be used to add empty lines which are ignored and help improve readability of the script).", ""]),
      ([RUN_COMMAND, "", "Save the script to file and then immediately run the script using the BASH command \"feasst < [script].txt > [script].log 2>&1\", where feasst is the compiled executable and [script].txt is the text file with the script. If the script has been previously-saved, then that file name will be used automatically. Otherwise, you will be prompted for the value of [script]."]),
    ]
    every_mc = copy.deepcopy(every_mc_no_for)
    #every_mc.insert(0, class_["For"])

    # Put RandomMT19937 first
    for iran,ran in enumerate(data['Random']):
        if ran[0] == 'RandomMT19937':
            data['Random'].insert(0, data['Random'].pop(iran))
            break

    data = data | {
    'tutorial': tutorial,
    'selected_index': 0,
    'scroll_des': 0,
    'node_key': UNIQUE_START,
    'script': "",
    'script_file': "",
    'num_arg_rows': 5,
    'num_script_rows': 5,
    'mc_stage': 0,
    'particles': particles,
    'particle_id': [],
    'particle_name': [],
    'configuration_name': [],
    'configuration_descript': [],
    'prepend': "",
    'group_name':[],
    'group_descript': [],
    'group_describing': False,
    'group_to_process': [],
    'excluded_args': [],
    'factory': factory,
    'previous_state': None,
    'future_state': None,
    'nested_options' : None,
    'user_message': '',
    UNIQUE_START: [
        ("[I'm ready to simulate now]", "simulate", """The selected option is shown highlighted at the top, and the description of your current selection is provided below.

Press the arrow keys (or w-a-s-d) to change your current selection.
Press the "enter" key to choose your current selection.
Press the "u" key to undo your last choice.
Press the "r" key to redo your last undo.

A mouse scroll wheel and swiping up, down, left or right on a touch screen or track pad may also function as arrow key presses.

Have you FEASST'ed with us before?

⠀⠀⣠⣾⣿⣿⣿⣷⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠀⣿⣿⣿⣿⣿⣿⣿⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣰⣶⣆⠀⠀⠀⠀
⠀⠀⢿⣿⣿⣿⣿⣿⡿⠀⠀⠀⠀⠀⠀⠀⠀⣴⣿⣿⣿⣿⣿⣦⠀⠀
⠀⠀⠀⠙⠛⠛⠛⠋⠀⠀⠀⠀⠀⠀⠀⢠⣤⣿⣿⣿⣿⣿⣿⣿⣤⡄
⢠⣤⣶⠀⠀⠀⠀⠀⣶⣶⣤⣤⡀⠀⠀⠸⠿⠿⠿⠿⠿⠿⠿⠿⠿⠇
⢸⣿⣿⡇⠀⠀⠀⢸⣿⣿⣿⣿⣷⡄⠀⠀⠀⣤⣤⣤⣤⣤⣤⠄⠀⠀
⢸⣿⣿⣿⡀⠀⠀⣾⣿⣿⣿⣿⣿⣿⣆⠀⣼⣿⣿⠛⠉⠉⠀⠀⠀⠀
⢸⣿⣿⣿⡇⠀⢰⣿⣿⣿⣿⡿⢿⣿⣿⣿⣿⡿⠁⠀⠀⠀⠀⠀⠀⠀
⢸⣿⣿⣿⣷⠀⢸⣿⣿⣿⣿⡇⠀⠹⣿⣿⡿⠁⠀⠀⠀⠀⠀⠀⠀⠀
⢸⣿⣿⣿⣿⣇⣿⣿⣿⣿⣿⡇⠀⠀⠈⠉⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀
⢸⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
⠀⠛⠿⠿⠿⠿⠿⠿⠿⠿⠛⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀

Like a New Jersey diner, there are a lot of FEASST options. Press the "enter" key to see the menu.

While I'm taking your order, I'm writing the text-based script needed to start a fresh Monte Carlo simulation. The simulation won't actually start cooking until that script is given to the executable (e.g., "feasst < script.txt" in BASH/etc).

You can also ask for things that aren't on the menu by editting a text-based script directly. And you can find a lot of FEASST recipes here: https://pages.nist.gov/feasst/tutorial/README.html , which also include post processing and analysis of the simulations."""),
        ("[What is FEASST?]", "", """The Free Energy and Advanced Sampling Simulation Toolkit (FEASST) is a free and open-source software to conduct molecular- and particle-based simulations with Monte Carlo methods.

New users can start with the website (https://pages.nist.gov/feasst/ or https://doi.org/10.18434/M3S095>), manuscript (https://doi.org/10.1063/5.0224283), GitHub discussion (https://github.com/usnistgov/feasst/discussions) and a five minute video (https://www.nist.gov/video/how-use-feasst-0255-monte-carlo-molecular-simulation-software).

Support FEASST with a GitHub (https://github.com/usnistgov/feasst) star or manuscript (https://doi.org/10.1063/5.0224283) citation!"""),
#        ("[What is this software?]", "", """This software is used to generate FEASST text scripts that instruct the FEASST executable on how to run the simulation. No simulation is actually run until the text script is given to the executable.
#
#This is the most beginner-friendly user interface (UI) available for the FEASST molecular simulation software. The drawbacks of an interactive UI, such as this one, compared to text-based scripts is reduced reproducibility and reduced automation. It can be hard to remember exactly which options were chosen in a UI and in which order, but a text file or script leaves a clear record of the recipe to exactly reproduce what you have done in the past. Similarly, text-based interfaces can be accessed programmatically to enable automation. Without automation, using a UI for many simulations becomes tedious.
#
#The major benefit of a UI is that it is easier to get started. While a UI is best for first time users, a text interface is better for intermediate to advanced users. The best of both worlds is to have a UI for new users that demonstrates how to move toward using the more advanced text interface. That is the goal of this UI. This UI describes each available option and builds the text script along the way. It enforces the order of initialization (e.g. Configuration, then Potential, then ThermoParams, etc) and it lists the only available options (e.g., selecting a particle for a move lists only those particles that were previously initialized).
#
#The main result of using this UI is a FEASST text script that the FEASST executable reads to run a simulation. After learning more about FEASST, editting text scripts directly instead of using this UI may be a better option. Thus, the ultimate goal of this UI is for the user to no longer need it."""),
        #("[Speak to a real person]", "", "Please contact the developers (see https://pages.nist.gov/feasst/CONTACT.html)"),
        #("[Who contributed to FEASST?]", "", "There is a growing community of FEASST users. Some are listed on the Acknowledgement page https://pages.nist.gov/feasst/dev/sphinx/ACKNOWLEDGEMENT.html"),
        #("[How can I contribute?]", "", "There are many ways to contribute, including contacting the developers with suggestions, issue or comments. You can also cite FEASST in your manuscripts, presentations, posters and other acknowledgement sections. Or you can contribute directly to the open-source code, share a new tutorial, plugin, model or other work helpful to FEASST users."),
        #("[How do I cite FEASST?]", "", "Citation of not only FEASST, but also the original papers of the specific methods that you used is recommended. See https://pages.nist.gov/feasst/dev/sphinx/CITATION.html for citing FEASST."),
        #("[Install FEASST]", "", "See https://pages.nist.gov/feasst for installation instructions."),
        (LOAD_SAVE_FILE, "", """Load a save file that was previously created with this program with a json extension. Type the [name] of an existing \"[name].json\" file without typing the json file extension. Type into the box and then press the ENTER key to submit."""),
        (LOAD_SCRIPT_FILE, "", """Load an existing script file. Type the name into the box and then press the ENTER key to submit."""),
        (OPEN_HTML, "", "Open the local HTML documentation specific to your exact version of FEASST, which is located on your computer in /feasst/data/html/index.html. May require \"sudo apt install xdg-utils\". For WSL users, \"sudo apt install wslu\" to enable xdg-open."),
        ("[Settings]", "settings", "Adjust the settings for this program."),
    ],
    "exit": [],
    "settings": [
        (NUMBER_ARG_ROWS, 'settings', "Set an integer number of rows that the selectable arguments at the top occupy.", ""),
        (NUMBER_SCRIPT_ROWS, 'settings', "Set an integer number of rows for the script preview window.", ""),
        (TOGGLE_TUTORIAL, 'settings', "Toggle tutorial mode on or off by inputting \"0\" (or \"False\") for off and \"1\" (or \"True\") for on.", ""),
        ("[Back]", UNIQUE_START, "Go back to the main menu.")
    ],
    "simulate": [
        ("MonteCarlo", "monte_carlo0", """Conduct a serial Monte Carlo simulation, post process or evaluate a reference Configuration.

This is the typical choice for most FEASST scripts."""),
        class_["Prefetch"],
        ("Server", "server", """Server is not yet documented here.

See https://pages.nist.gov/feasst/plugin/server/README.html"""),
        ("Restart", "Restart", """Restart a simulation from a previously-created Checkpoint file."""),
    ],

# mc_stage: stages of MonteCarlo
# 0: random not set
# 1: config not set
# 2: potential not set
# 3: thermo param not set
# 4: criteria not set
# 5: criteria set
    "monte_carlo0": [
        #class_["RandomMT19937"],
        #class_["RandomModulo"],
        class_["Random"],
        class_config,
    ] + every_mc,
    "monte_carlo1": [
        class_config,
    ] + every_mc,
    "monte_carlo2": [
        class_["Potential"],
        class_["RefPotential"],
        class_config,
    ] + every_mc,
    "monte_carlo3": [
        class_['ThermoParams'],
        class_["Potential"],
        class_["RefPotential"],
    ] + every_mc,
    "monte_carlo4": [
        class_['Criteria'],
        class_['ThermoParams'],
    ] + every_mc_no_for,
    "monte_carlo5": [
        class_['Trial'],
        class_['Analyze'],
        class_['Modify'],
        class_['Action'],
        class_['Criteria'],
        class_['ThermoParams'],
    ] + every_mc,
}
    # alphabetize data keys so that changes on different computers are less likely to bloat revision control
    data = dict(sorted(data.items()))
    return data

def enter_to_submit_validator(ch):
    if ch in (10, curses.KEY_ENTER):  # Enter key
        return 7  # Ctrl-G code to signal end of editing
    return ch

def str_to_bool(value: str) -> bool:
    """
    Convert a string to a boolean.
    Accepts common truthy/falsey values (case-insensitive).
    Raises ValueError for unrecognized inputs.
    """
    if not isinstance(value, str):
        raise TypeError("Input must be a string")

    truthy = {"True", "true", "yes", "y", "1", "on"}
    falsey = {"False", "false", "no", "n", "0", "off"}

    val = value.strip().lower()
    if val in truthy:
        return True
    elif val in falsey:
        return False
    else:
        raise ValueError(f"Invalid boolean string: '{value}'")

def input_box(stdscr, width, instruct, data, default_text=''):
# TESTING     if stdscr == None:
# TESTING         return ""
# TESTING     if isinstance(stdscr, list):
# TESTING         if len(stdscr) == 0:
# TESTING             return ""
# TESTING         return stdscr.pop(0)

    curses.curs_set(2) # turn on cursor
    stdscr.addstr(0, 0, instruct+str((width-len(instruct))*' '))

    edit_height, edit_width = data['num_arg_rows']-2, width-2
    edit_y, edit_x = 2, 1
    editwin = curses.newwin(edit_height, edit_width, edit_y, edit_x)

    # Draw a rectangle border around the textbox
    rectangle(stdscr, edit_y - 1, edit_x - 1,
              edit_y + edit_height, edit_x + edit_width)

    stdscr.refresh()

    editwin.addstr(0, 0, default_text)
    editwin.move(0, len(default_text))
    editwin.refresh()

    box = Textbox(editwin)
    box.edit(enter_to_submit_validator)
    text = box.gather().strip().replace('\n','')
    curses.curs_set(0) # hide the blinking cursor
    return text

def print_bar(width, label=''):
    halfwidth = int(width/2-len(label)/2)
    return str(halfwidth*'=')+label+str((width-halfwidth-len(label))*'=')

def update_mc_stage(data, chc0, stdscr):
    if chc0[:6] == 'Random':
        assert data['mc_stage'] == 0, data['mc_stage']
        data['mc_stage'] = 1
    elif chc0 == 'Configuration':
        assert data['mc_stage'] <= 2, data['mc_stage']
        data['mc_stage'] = 2
    elif chc0 == 'Potential':
        assert data['mc_stage'] == 2, data['mc_stage']
        data['mc_stage'] = 3
    elif chc0 == 'ThermoParams':
        assert data['mc_stage'] >= 3, data['mc_stage']
        if data['mc_stage'] == 3:
            data['mc_stage'] = 4
    elif any(crit == chc0 for crit in data['factory']['Criteria']):# == 'Metropolis' or chc0 == 'MayerSampling' or chc0 == 'FlatHistogram':
        assert data['mc_stage'] >= 4, data['mc_stage']
        data['mc_stage'] = 5

def wrap_in_width(text, width, join):
    wrapped_lines = []
    for line in text.splitlines():
        if line.strip():  # Non-empty line
            wrapped_lines.extend(textwrap.wrap(line, width=width-1)) # avoid a line that is just one space
        else:
            wrapped_lines.append("")  # Preserve blank lines
    if join:
        return "\n".join(wrapped_lines)
    else:
        return wrapped_lines

def user_message(message, stdscr, title='WARNING'):
    height, width = stdscr.getmaxyx()
    box_height = height - 2
    box_width = width - 2
    start_y = (height - box_height) // 2
    start_x = (width - box_width) // 2
    win = curses.newwin(box_height, box_width, start_y, start_x)
    win.box()
    win.addstr(0, (box_width - len(title)) // 2, title, curses.A_BOLD | curses.A_REVERSE)
    border=2
    #win.addstr(2, 2, message, curses.A_BOLD)
    for iline,line in enumerate(wrap_in_width(message, width=box_width-border-1, join=False)):
        if iline < box_height - 6:
            win.addstr(border+iline, border, line, curses.A_BOLD)
    if title == "TUTORIAL":
        instruction = "Disable tutorials in main menu [Settings]"
        win.addstr(box_height - 4, (box_width - len(instruction)) // 2, instruction, curses.A_DIM)
        instruction = "or Python argument \"--disable-tutorial\""
        win.addstr(box_height - 3, (box_width - len(instruction)) // 2, instruction, curses.A_DIM)
    instruction = "Press any key to continue..."
    win.addstr(box_height - 2, (box_width - len(instruction)) // 2, instruction, curses.A_BOLD | curses.A_REVERSE)
    win.getch()

def tutorial(data, stdscr, message):
    if data['tutorial']:
        user_message(message, stdscr, title='TUTORIAL')

def display(data, stdscr, choices):
    log.debug('display.node_key:'+str(data['node_key']))
    height, width = stdscr.getmaxyx()
    # first, display the arguments with data['num_arg_rows']
    # the number of columns allowed depends on the width
    selected_index = data['selected_index']
    sel = selected_index
    num_rows = data['num_arg_rows']
    assert type(num_rows) == int
    if num_rows <= 2:
        user_message("The number of rows must be greater than 2. Setting to the minimum of 3.", stdscr=stdscr)
        data['num_arg_rows'] = 3
    min_col_width = 3+max(len(row[0]) for row in choices)
    log.debug('min_col_width'+str(min_col_width))
    num_cols = int(width/min_col_width)
    col_width = int(width/num_cols)

    # display arguments
    if len(choices) <= num_rows*num_cols:
        chcs = choices
    else:
        half_row_up = int((num_rows+1)/2)
        half_row_down = int((num_rows-1)/2)
        half_col_up = int((num_cols+1)/2)
        half_col_down = int((num_cols-1)/2)
        half_col = half_col_down
        if num_cols % 2 == 0:
            half_col = half_col_up
        if selected_index < half_col_down*num_rows + half_row_down:
            chcs = choices[:num_rows*num_cols]
        elif selected_index > len(choices) - half_col*num_rows - half_row_up:
            chcs = choices[len(choices)-num_rows*num_cols:]
            sel = num_rows*num_cols-(len(choices) - selected_index)
        else:
            chcs = choices[(selected_index-half_row_down-half_col_down*num_rows):
                           (selected_index+half_row_down*half_col_up*num_rows)]
            sel = half_col_down*num_rows + half_row_up - 1
    for ir in range(num_rows):
        for ic in range(num_cols):
            idx=ic*num_rows+ir
            if idx < len(chcs):
                element = chcs[idx]
                choice_text = element[0]
                extra=col_width - 3 - len(choice_text)
                if ic == num_cols - 1:
                    extra -= 1
                if idx == sel:
                    stdscr.addstr(f"> {choice_text} "+str(extra*' '), curses.color_pair(2))
                else:
                    stdscr.addstr(f"  {choice_text} "+str(extra*' '), curses.color_pair(1))
        stdscr.addstr('\n')

    # display script, if exists
    lines_for_descript = height - (data['num_arg_rows'] + 1)
    if data['script'] != '':
        stdscr.addstr(print_bar(width=width, label='SCRIPT'))
        lines = data['script'].split('\n')
        if lines[-1] == '':
            lines.pop()
        text = ''
        if len(lines) > data['num_script_rows']:
            text += '...\n'
            for line in lines[-data['num_script_rows']+1:]:
                text += line+'\n'
        else:
            for line in lines[-data['num_script_rows']:]:
                text += line+'\n'
        stdscr.addstr(text)
        lines_for_descript -= min(len(lines), data['num_script_rows'])+1

    # display description
    stdscr.addstr(print_bar(width=width, label='DESCRIPTION'))
    text = choices[selected_index][2]
    # wrap text while perserving endlines
    wrapped_lines = wrap_in_width(text=text, width=width, join=False)
    if len(wrapped_lines) <= lines_for_descript:
        stdscr.addstr("\n".join(wrapped_lines))
    else:
        scroll_des = data['scroll_des']
        if scroll_des > 0:
            stdscr.addstr('[Press "q" to scroll up]\n', curses.color_pair(2))
            lines_for_descript -= 1
        display_down_instruct = False
        if len(wrapped_lines) - scroll_des > lines_for_descript:
            stdscr.addstr("\n".join(wrapped_lines[scroll_des:scroll_des+lines_for_descript-1]))
            stdscr.addstr('\n[Press "e" to scroll down]\n', curses.color_pair(2))
        else:
            stdscr.addstr("\n".join(wrapped_lines[scroll_des:]))
            data['scroll_des'] = len(wrapped_lines) - lines_for_descript -1

    if data['user_message'] != '':
        user_message(message=data['user_message'], stdscr=stdscr)
        data['user_message'] = ''

def write_files(data):
    with open(data['script_file']+'.txt', 'w', encoding='utf-8') as file1:
        file1.write(data['script'])
    data['previous_state'] = None
    data['future_state'] = None
    with open(data['script_file']+'.json', 'w', encoding='utf-8') as file1:
        json.dump(data, file1, indent=2)

def update(choices, data, stdscr):
    selected_index = data['selected_index']
    width = None
    if stdscr != None:
        height, width = stdscr.getmaxyx()
    def node_key(): return data['node_key']

    chc = choices[selected_index]
    log.debug("chc0:" + chc[0] + " chc1:" + chc[1] + " selected_index:"+str(selected_index) + " len:" + str(len(data[node_key()])))

    # tutorial
    if chc[0] == 'Configuration':
        tutorial(data, stdscr, """In the next step, select the types of particles that may exist in the Configuration. Each particle type is described by a FEASST particle file. Even if the Configuration is to have 500 particles of the same type, that type of particle need only be named once here.

A particle type file should have only one particle each. If the Configuration is to have multiple types of particles, select separate particle type files one after the other.""")
    elif chc[1] == 'particles_select1':
        tutorial(data, stdscr, """In the next step, enter a name for the particle type. This name can later be used to refer to the particle types.

If no name is provided, the default name of the integer index of the particle type starts from \"0\".

Note that the path to the particle file in the resulting script begins with /feasst, which points to the feasst installation directory.""")
    if chc[0] == UNIQUE_PARTICLE:
        text = input_box(stdscr=stdscr, width=width, instruct="Type the name of the particle (Press Enter to submit):", data=data)
        data['node_key'] = chc[1]
    elif chc[1] in data:
        log.debug('chc1:'+str(chc[1])+' is in data')
        if len(chc) > 3:
            log.debug('parsing argument, not a class, and a value may need to be input')
            exclude_arg = True
            arg = chc[0]
            if chc[0] != NUMBER_ARG_ROWS and chc[0] != NUMBER_SCRIPT_ROWS and chc[0] != TOGGLE_TUTORIAL and chc[0] != COMMENT_LINE and chc[0] != CUSTOM_LINE and chc[0] != LOAD_PARTICLE_FILE and arg.find('[') != -1:
                log.debug('this argument has additional variables defined within [square brackets]')
                while arg.find('[') != -1:
                    param = arg[arg.find('[')+1 : arg.find(']')]
                    text = input_box(stdscr=stdscr, width=width, instruct="Type the value of ["+param+"] (Press Enter to submit):", data=data)
                    arg = arg.replace('['+param+']', text)
                    exclude_arg = False
            if chc[0] == LOAD_PARTICLE_FILE:
                path = input_box(stdscr=stdscr, width=width, instruct="Type the full path and file name of the particle file:", data=data)
                name = input_box(stdscr=stdscr, width=width, instruct="Type the name of the particle (Press Enter to submit):", data=data)
                if data['node_key'] == 'particles_select0':
                    data['script'] += ' '+data['prepend']+'particle_type='
                else:
                    data['script'] += ','
                data['script'] += name + ':' + path
                exclude_arg = False
                data['node_key'] = chc[1]
            elif data['node_key'][:-1] == 'particles_select':
                log.debug('obtain name and add to list of particle names')
                text = input_box(stdscr=stdscr, width=width, instruct="Type the name of the particle (Press Enter to submit):", data=data)
                if text == '':
                    text = str(len(data['particle_name']))
                if text in data['particle_name']:
                    pid = data['particle_name'].index(text)
                else:
                    data['particle_name'].append(text)
                    data['particle_id'].append(chc[0])
                    pid = -1
                if data['node_key'] == 'particles_select0':
                    data['script'] += ' '+data['prepend']+'particle_type='
                else:
                    data['script'] += ','
                data['script'] += data['particle_name'][pid] + ':/feasst/' + chc[3]
                log.debug('node:'+str(data['node_key']))
                exclude_arg = False
                data['node_key'] = chc[1]
            elif chc[3] == "factory":
                data['nested_options'] = data['node_key']
                data['node_key'] = chc[1]
                data['script'] += " "+chc[0]+"="
            elif chc[3] == "particle_type_chooser":
                assert len(data['particle_name']) > 0, "No particles defined in Configuration."
                assert len(data['particle_name']) == len(data['particle_id'])
                data['particle_type_chooser'] = []
                for idx,part in enumerate(data['particle_name']):
                    pid = data['particle_id'][idx]
                    descript = data['particles'][pid][0]
                    data['particle_type_chooser'].append((part, chc[1], descript))
                data['script'] += ' '+data['prepend']+'particle_type='
                data['node_key'] = 'particle_type_chooser'
            elif chc[3] == "configuration_chooser":
                if len(data['configuration_name']) > 0: return#, "No configuration names defined."
                assert len(data['configuration_name']) == len(data['configuration_descript'])
                data['configuration_chooser'] = []
                for idx,config in enumerate(data['configuration_name']):
                    data['configuration_chooser'].append((config, chc[1], data['configuration_descript'][idx]))
                data['script'] += " config="
                data['node_key'] = 'configuration_chooser'
            elif chc[3] == "group_chooser":
                if len(data['group_name']) == 0:
                    user_message('No groups were defined in Configuration.', stdscr)
                    return
                data['group_chooser'] = []
                assert len(data['group_name']) == len(data['group_descript']), str(data['group_name'])+' '+str(data['group_descript'])
                for idx,group in enumerate(data['group_name']):
                    data['group_chooser'].append((group, chc[1], data['group_descript'][idx]))
                data['script'] += " group="
                data['node_key'] = 'group_chooser'
            else:
                log.debug('parsing argument')
                text = input_box(stdscr=stdscr, width=width, instruct="Type the value of \""+arg+"\" (Press Enter to submit):", data=data)
                if chc[0] == NUMBER_ARG_ROWS:
                    exclude_arg = False
                    try:
                        data['num_arg_rows'] = int(text)
                    except:
                        user_message("The number of rows must be an integer.", stdscr)
                elif chc[0] == NUMBER_SCRIPT_ROWS:
                    exclude_arg = False
                    try:
                        data['num_script_rows'] = int(text)
                    except:
                        user_message("The number of rows must be an integer.", stdscr)
                elif chc[0] == TOGGLE_TUTORIAL:
                    exclude_arg = False
                    truthy = {"true", "t", "yes", "y", "1", "on"}
                    falsey = {"false", "f", "no", "n", "0", "off"}
                    if text in truthy:
                        data['tutorial'] = True
                    elif text in falsey:
                        data['tutorial'] = False
                    else:
                        user_message("The tutorial mode must be entered as 0, False, 1 or True.", stdscr)
                    if data['tutorial']:
                        user_message('The tutorial mode is ON', stdscr=stdscr, title='TUTORIAL')
                else:
                    data['script'] += " " + data['prepend'] + arg + "=" + text
                    if chc[3] != "" and chc[3] in data:
                        data['node_key'] = chc[3]
                log.debug('here1234:'+chc[0]+' '+chc[1])
                if data['group_describing']:
                    data['group_descript'][-1] += arg + "=" + text
                if chc[0] == 'name' and chc[1] == 'Configuration'+OPTIONAL:
                    data['configuration_name'].append(text)
                elif chc[0] == 'group' and data['node_key'] == 'Configuration'+OPTIONAL:
                    data['nested_options'] = data['node_key']
                    log.debug(data['node_key']+'chc1:'+chc[1])
                    data['node_key'] = chc[1]
                    grps = text.split(',')
                    data['prepend'] = grps[0]+'_'
                    data['group_name'].append(grps[0])
                    for grp in grps[1:]:
                        if grp != "":
                            data['group_to_process'].append(grp)
                    data['user_message'] = "Now taking arguments for group \""+data['prepend'][:-1]+"\""
                    data['group_descript'].append("")
                    data['group_describing'] = True
            # remove from choices (unless subsitutable variable)
            if exclude_arg:
                data['excluded_args'].append(arg)
        else:
            # this is a class or special option, not an argument
            log.debug("else")
            data['node_key'] = chc[1]
            if chc[0] == "MonteCarlo":
                data['script'] += "MonteCarlo\n"
                tutorial(data, stdscr, """A MonteCarlo simulation is initialized in a particular order. The first class is typically either a Random number generator, or a Configuration. This is the only stage that a Random number generator can be chosen, or else the default RandomMT19937 is used with the current time as the seed.""")#+str(10000*'-'))
            elif chc[0] == "Restart":
                data['script'] += "Restart"
                data['mc_stage'] = 5
            elif chc[0] == "Prefetch":
                data['script'] += "Prefetch"
            else:
                if data['script'] != "":
                    if chc[0][:len(NO_MORE_ARGUMENTS)] == NO_MORE_ARGUMENTS:
                        log.debug('parsing no more arguments')
                        data['prepend'] = ""
                        if data['node_key'] == 'Configuration'+OPTIONAL:
                            log.debug('record config line')
                            lastline = data['script'].split('\n')[-1]
                            data['configuration_descript'].append("Below is the line that created this Configuration name (up until the name was entered):\n\n"+lastline)
                            if len(data['configuration_name']) == len(data['configuration_descript']) - 1:
                                data['configuration_name'].append(str(len(data['configuration_name'])))
                            if data['mc_stage'] == 2:
                                tutorial(data, stdscr, "Now that a Configuration is initialized, another Configuration could be created for the Gibbs ensemble, or you can move on to describing the Potential energy functions for interactions between particles. This is the last chance to add another Configuration.")
                        if data['node_key'] == 'Potential'+OPTIONAL:
                            tutorial(data, stdscr, "Now that a Potential is initialized, another Potential could be created, or you can move on to describing the ThermoParams. This is the last chance to add another Potential.")
                        if data['node_key'] == 'ThermoParams'+OPTIONAL:
                            tutorial(data, stdscr, "Now that ThermoParams is initialized, it is time to select what kind of acceptance Criteria you will use in your MonteCarlo simulation. You can change this acceptance Criteria by making a new one at any time.")
                        if any(crit+OPTIONAL == data['node_key'] for crit in data['factory']['Criteria']):
                            tutorial(data, stdscr, "Now that the acceptance Criteria is initialized, all Trial, Action, Analyze, Modify and Criteria are available as options.")
                        if len(data['group_to_process']) > 0:
                            data['prepend'] = data['group_to_process'][0] + '_'
                            data['group_name'].append(data['group_to_process'].pop(0))
                            data['user_message'] = "Now taking arguments for group \""+data['prepend'][:-1]+"\""
                            data['excluded_args'] = []
                            data['group_descript'].append("")
                            data['group_describing'] = True
                        else:
                            data['group_describing'] = False
                            if data['nested_options']:
                                data['node_key'] = data['nested_options']
                                data['nested_options'] = None
                            else:
                                data['script'] += '\n'
                                data['node_key'] = 'monte_carlo'+str(data['mc_stage'])
                                data['excluded_args'] = []
                    elif chc[0] == NO_MORE_PARTICLES:
                        log.debug('parsing no more particles')
                        tutorial(data, stdscr, "In the next step, the Domain must be defined by one of the following arguments.")
                    else:
                        if chc[0] in data['factory']:
                            log.debug('this is a factory, no need to print or update stage')
                        else:
                            log.debug('this is a class')
                            data['script'] += chc[0]
                            update_mc_stage(data, chc[0], stdscr)
                            if chc[0] == 'Checkpoint':
                                data.update(remove_checkpoint(data))

    elif chc[0] == LOAD_SAVE_FILE:
        log.debug("Loading save file")
        instruct = "Type file name (Press ENTER key to submit):"
        flname = input_box(stdscr=stdscr, width=width, instruct=instruct, data=data) + '.json'
        if os.path.exists(flname):
            with open(flname, 'r') as file1:
                data.update(json.load(file1))
        else:
            user_message("The file:\""+flname+"\" does not exist. Try again.", stdscr)
    elif chc[0] == LOAD_SCRIPT_FILE:
        log.debug('Loading script file')
        instruct = "Type file name (Press ENTER key to submit):"
        flname = input_box(stdscr=stdscr, width=width, instruct=instruct, data=data)
        if os.path.exists(flname):
            with open(flname, 'r') as file1:
                dflt_data = str(files('feasst').joinpath('data/menu.json'))
                data.update(upload_script(file1, dflt_data))
            #    data.update(json.load(file1))
        else:
            user_message("The file:\""+flname+"\" does not exist. Try again.", stdscr)
    elif chc[0] == OPEN_HTML:
        file_path = str(files('feasst').joinpath('data/html/index.html'))
        webbrowser.open(f"file://{file_path}")
    elif chc[0] == SAVE_FILE:
        log.debug("Saving")
        tutorial(data, stdscr, "Saving creates two files. The first is the script with a .txt file extension. The second is the current state of the menu with a .json file extension, which can be loaded with the command \"feasst-menu --load file_name.json\". Enter the name of both of these files without including the file extension.")
        instruct = "Type file name (Press ENTER key to submit):"
        default_text = data['script_file']
        data['script_file'] = input_box(stdscr=stdscr, width=width, instruct=instruct, data=data, default_text=default_text)
        if data['script_file'] != "":
            write_files(data)
        log.debug("script_file:" +data['script_file'])
    elif chc[0] == COMMENT_LINE:
        text = input_box(stdscr=stdscr, width=width, instruct="Type the script comment (Press Enter to submit):", data=data)
        exclude_arg = False
        data['script'] += '# '+text+'\n'
        data['node_key'] = 'monte_carlo'+str(data['mc_stage'])
    elif chc[0] == CUSTOM_LINE:
        text = input_box(stdscr=stdscr, width=width, instruct="Type the custom line in script (Press Enter to submit):", data=data)
        exclude_arg = False
        data['script'] += text+'\n'
        data['node_key'] = 'monte_carlo'+str(data['mc_stage'])
    elif chc[0] == RUN_COMMAND:
        data['script_file'] = input_box(stdscr=stdscr, width=width, instruct="Type file name (Press ENTER key to submit)", data=data, default_text=data['script_file'])
        write_files(data)
        syscode = subprocess.call("feasst < "+data['script_file']+'.txt > '+data['script_file']+'.log 2>&1', shell=True, executable='/bin/bash')
        if syscode != 0:
            user_message("An error was encountered while running the script. Please see recent file output for the error message.", stdscr)
        stdscr.refresh()
    elif selected_index < len(data[node_key()]):
        log.debug("Remove from choices")
        data[node_key()].pop(selected_index)

def update_previous_state(data):
    log.debug('store previous state and remove after so many previous states')
    data['previous_state'] = copy.deepcopy(data)
    # only enable so many previous states to be saved to reduce memory requirements.
    ps = 'previous_state'
    if data[ps][ps]:
        if data[ps][ps][ps]:
            if data[ps][ps][ps][ps]:
                if data[ps][ps][ps][ps][ps]:
                    if data[ps][ps][ps][ps][ps][ps]:
                        if data[ps][ps][ps][ps][ps][ps][ps]:
                            if data[ps][ps][ps][ps][ps][ps][ps][ps]:
                                if data[ps][ps][ps][ps][ps][ps][ps][ps][ps]:
                                    data[ps][ps][ps][ps][ps][ps][ps][ps][ps] = None
#                                    if data[ps][ps][ps][ps][ps][ps][ps][ps][ps][ps]:
#                                        if data[ps][ps][ps][ps][ps][ps][ps][ps][ps][ps][ps]:
#                                            if data[ps][ps][ps][ps][ps][ps][ps][ps][ps][ps][ps][ps]:
#                                                if data[ps][ps][ps][ps][ps][ps][ps][ps][ps][ps][ps][ps][ps]:
#                                                    if data[ps][ps][ps][ps][ps][ps][ps][ps][ps][ps][ps][ps][ps][ps]:
#                                                        if data[ps][ps][ps][ps][ps][ps][ps][ps][ps][ps][ps][ps][ps][ps][ps]:
#                                                            data[ps][ps][ps][ps][ps][ps][ps][ps][ps][ps][ps][ps][ps][ps][ps] = None

# TESTING def run(stdscr, data, repro=None):
def run(stdscr, data):
# TESTING     testing = False
# TESTING     if isinstance(repro, list):
# TESTING         testing = True
# TESTING         print('parsing', repro)
# TESTING         if len(repro) == 0:
# TESTING             return
# TESTING     print('testing', testing, type(repro))

# TESTING     if not testing:
    curses.curs_set(0) # hide the blinking cursor

    # Initialize color pairs (White text on black background, and highlighted text)
    curses.init_pair(1, curses.COLOR_WHITE, curses.COLOR_BLACK)
    curses.init_pair(2, curses.COLOR_BLACK, curses.COLOR_WHITE)

    def node_key(): return data['node_key']
#    stdscr.erase()
#    try:
#        if data['tutorial']:
#            #user_message("hi", stdscr=stdscr)
#            user_message(data[UNIQUE_START][0][2], stdscr=stdscr)
#    except curses.error:
#        pass

    while True:
# TESTING         testing = False
# TESTING         if isinstance(repro, list):
# TESTING             testing = True
# TESTING             print('parsing', repro)
# TESTING             if len(repro) == 0:
# TESTING                 data['script_file'] = 'test'
# TESTING                 write_files(data)
# TESTING                 return
# TESTING         print('testing', testing, type(repro))
# TESTING
# TESTING         if not testing:
        stdscr.erase()
        height, width = stdscr.getmaxyx()

        # handle exitting
        if node_key() == "exit":
            return

        choices = []
        num_while = 0
        while len(choices) <= 1:
        #while len(data[node_key()]) <= 1:
            choices = copy.deepcopy(data[node_key()])
            log.debug('excluded args: '+str(data['excluded_args']))
            #log.debug('Remove excluded args from '+str(str(choices)))
            # remove excluded args
            for idx in reversed(range(len(choices))): #,chc in enumerate(choices):
                chc=choices[idx]
                log.debug('chc0:'+str(chc[0]))
                if chc[0] in data['excluded_args']:
                    choices.pop(idx)
            #log.debug('Removed excluded args from '+str(str(choices)))

            num_while += 1
            if num_while > 1e3:
                user_message("An infinite loop was obtained. Please report how to reproduce this behavior to the developers.", stdscr=stdscr)
                return
            if len(choices) <= 1:
                log.debug('automatically select choice if only one or zero choices')
            if len(choices) == 1:
                log.debug('selecting only option')
                data['selected_index'] = 0
                display(data=data, stdscr=stdscr, choices=choices)
                update(choices=choices, data=data, stdscr=stdscr)
                data['scroll_des'] = 0
                stdscr.erase()
            elif len(choices) == 0:
                log.debug(node_key()+OPTIONAL)
                if node_key()+OPTIONAL in data:
                    log.debug('moving to optional')
                    #update(selected_index=-1, choices=data['node_key'], data=data, stdscr=stdscr)
                    data['node_key'] = node_key() + OPTIONAL
                    data['selected_index'] = 0 # reset selection for the new node
                    data['scroll_des'] = 0
                else:
                    assert False
        #if num_while > 1:
        #    continue

        # add default choices
        defaults = []
        if data['script'] != "":
            defaults += [(SAVE_FILE, "", data['script'])]
        defaults += [("[Exit]", "exit", "Exit the program.")]
        choices += defaults

        # render text
        try:
# TESTING             if not testing:
            display(data=data, stdscr=stdscr, choices=choices)
        except curses.error:
            pass

# TESTING         if testing:
# TESTING             print('repro', repro)
# TESTING             if len(repro):
# TESTING                 return
# TESTING             data['selected_index'] = int(repro.pop(0))
# TESTING             update(choices=choices, data=data, stdscr=None)
# TESTING         else:
        stdscr.refresh()

        # capture keyboard input
        key = stdscr.getch()
        if key == curses.KEY_UP or key == ord('w'):
            data['selected_index'] = (data['selected_index'] - 1) % len(choices)
        elif key == curses.KEY_DOWN or key == ord('s'):
            data['selected_index'] = (data['selected_index'] + 1) % len(choices)
        elif key == curses.KEY_LEFT or key == ord('a'):
            data['selected_index'] -= data['num_arg_rows']
            if data['selected_index'] < 0:
                data['selected_index'] = 0
        elif key == curses.KEY_RIGHT or key == ord('d'):
            data['selected_index'] += data['num_arg_rows']
            if data['selected_index'] >= len(choices):
                data['selected_index'] = len(choices) - 1
        elif key == ord('q'):
            if data['scroll_des'] > 0:
                data['scroll_des'] -= 1
        elif key == ord('e'):
            data['scroll_des'] += 1
        elif key == ord('u'):
            if data['previous_state']:
                old_data = copy.deepcopy(data)
                data = copy.deepcopy(data['previous_state'])
                data['future_state'] = old_data
            else:
                tutorial(data, stdscr, "There are no more previous states saved for further undo operations. Only a finite number of previous states are saved to reduce memory requirements.")
        elif key == ord('r'):
            if data['future_state']:
                data = copy.deepcopy(data['future_state'])
            else:
                tutorial(data, stdscr, "There are no more future states saved for further redo operations. Only a finite number of future states are saved to reduce memory requirements. In addition, future states are cleared upon pressing the enter key.")
        elif key == ord('\n'): # ENTER key
            update_previous_state(data)
            data['future_state'] = None
            update(choices=choices, data=data, stdscr=stdscr)
            data['selected_index'] = 0 # reset selection for the new node
            data['scroll_des'] = 0

def get_mc_stage(data):
    script = data['script']
    stage = 0
    if 'Random' in script:
        stage = 1
    if 'Configuration' in script:
        stage = 2
    if 'Potential' in script:
        stage = 3
    if 'ThermoParams' in script:
        stage = 4
    if 'Metropolis' in script or 'MayerSampling' in script or 'FlatHistogram' in script:
        stage = 5
    return stage

def get_particles(data):
    ids = list()
    names = list()
    for line in data['script'].split('\n'):
        #print('line:', line, '|', line[:13])
        if line != '':
            if line[:13] == 'Configuration':
                #print('here:', line)
                for arg in line.split(' '):
                    #print('arg', arg)
                    pairs=arg.split('=')
                    #print('pairs:', pairs)
                    if pairs[0] == 'particle_type':
                        #print('types:', pairs[1])
                        for typ in pairs[1].split(','):
                            #print('typ', typ)
                            p2 = typ.split(':')
                            assert len(p2) == 2
                            names.append(p2[0])
                            iid=p2[1].split('/')[-1]
                            ids.append(iid)
                            #print('iid:', iid)
                            #quit()
    #print(data['particles'])
    #print(data['particles'][ids[0]])
    #quit()
    return ids, names

def remove_checkpoint(data):
    index = 0
    while 'monte_carlo'+str(index) in data:
        node='monte_carlo'+str(index)
        for iopt,opt in reversed(list(enumerate(data[node]))):
            if data[node][iopt][0] == 'Checkpoint':
                data[node].pop(iopt)
        index += 1
        assert index < 1e3, "infinite loop"
    return data

def upload_script(file1, dflt_data):
    with open(dflt_data, 'r') as file2:
        data = json.load(file2)
    data['script'] = file1.read()
    data['mc_stage'] = get_mc_stage(data)
    data['node_key'] = 'monte_carlo'+str(data['mc_stage'])
    data['particle_id'], data['particle_name'] = get_particles(data)
    for line in data['script'].split('\n'):
        if line[:len(str('Checkpoint'))] == 'Checkpoint':
            data.update(remove_checkpoint(data))
    return data

def main_function():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--disable-tutorial", action="store_false", dest="tutorial", help="Disable tutorial")
    parser.add_argument('--load', '-l', type=str, default='', help='Optional json save file or script text file to load.')
    parser.add_argument('--descripts', '-d', type=str, default='', help='For developers only.')
    parser.add_argument('--particles', '-p', type=str, default='', help='For developers only.')
    parser.add_argument('--factory', '-f', type=str, default='', help='For developers only.')
    #parser.add_argument("--test_reproduction", type=str, default='', help='For developers only.')
    parser.set_defaults(tutorial=True) # Default is True
    args, unknown_args = parser.parse_known_args()
    assert len(unknown_args) == 0, 'An unknown argument was included: '+str(unknown_args)
    params = vars(args)
    particles = []

    # load data from file
    dflt_data = str(files('feasst').joinpath('data/menu.json'))
    if params['descripts'] == '' and params['particles'] == '' and params['factory'] == '':
        if params['load'] == '':
            params['load'] = dflt_data
        if os.path.exists(params['load']):
            with open(params['load'], 'r') as file1:
                if params['load'][-5:] == '.json':
                    data = json.load(file1)
                elif params['load'][-4:] == '.txt':
                    data = upload_script(file1, dflt_data)
                else:
                    assert False, 'unrecognized file extension in '+params['load']
        else:
            assert False, str(params['load']) + ' does not exist.'
# TESTING         if params['test_reproduction'] != '':
# TESTING             repro = params['test_reproduction'].split(',')
# TESTING             print(repro)
# TESTING             data['tutorial'] = False
# TESTING             run(None, data, repro)
# TESTING         else:
        data['tutorial'] = params['tutorial']
        curses.wrapper(run, data)
    else:
        # create the file that can be loaded
        # load particles
        if os.path.exists(params['particles']):
            with open(params['particles'], 'r') as file1:
                particles = json.load(file1)
        else:
            assert False, params['particles'] + ' does not exist.'

        # load descriptions
        if os.path.exists(params['descripts']):
            with open(params['descripts'], 'r') as file1:
                descript = json.load(file1)
        else:
            assert False, params['descripts'] + ' does not exist.'

        # load factory
        if os.path.exists(params['factory']):
            with open(params['factory'], 'r') as file1:
                factory = json.load(file1)
        else:
            assert False, params['factory'] + ' does not exist.'

        data = gen_data(descript=descript, particles=particles, factory=factory, tutorial=params['tutorial'])

        with open(params['load'], 'w', encoding='utf-8') as file1:
            json.dump(data, file1, indent=2)

if __name__ == "__main__":
    main_function()

