//#include <iostream>
#include "utils/include/debug.h"
#include "utils/include/io.h"
#include "utils/include/arguments.h"
#include "utils/include/arguments_extra.h"

namespace feasst {

std::string str(const std::string& key, const argtype& args) {
  auto pair = args.find(key);
  if (pair != args.end()) {
    return pair->second;
  }
  return std::string();
}

std::pair<std::string, argtype> parse_line(const std::string line,
  argtype * variables,
  bool * assign_to_list) {
  // variable replacement if not Let
  std::string new_line = line;
  if (variables) {
    for (const auto& pair : *variables) {
      std::stringstream ss3(line);
      std::string major3;
      ss3 >> major3;
      if (major3 != "Let") {
        replace(pair.first, pair.second, &new_line);
      }
    }
  }
  std::stringstream ss(new_line);
  std::string major;
  ss >> major;
  argtype args;
  if (major == "Let") {
    DEBUG("major:" << major);
    major = "";
    ASSERT(ss.str().find('=') != std::string::npos,
      "Let requires an \"=\" to define a variable value");
    ASSERT(ss.str().find('[') != std::string::npos,
      "Let requires a [ to enclose a [variable].");
    ASSERT(ss.str().find("]=") != std::string::npos,
      "Let requires a \"]=\" to enclose a [variable] and define its value.");
    std::stringstream ss2(new_line);
    std::string variable;
    std::getline(ss2, variable, '[');
    std::getline(ss2, variable, ']');
    //std::cout << "variable:" << variable << std::endl;
    std::string value;
    std::getline(ss2, value, '=');
    std::getline(ss2, value);
    //std::cout << "value:" << value << std::endl;
    //DEBUG("value:"<<value);
    auto pair = variables->find(variable);
    if (pair != variables->end()) {
      // Allow identical Let for monte_carlo tutorial 2 (multi-mc)
      ASSERT(pair->second == value, "variable:\"" << variable << "\" is already defined with a different value:" << value);
    }
    (*variables)["["+variable+"]"] = value;
  } else {
    int num_pairs = 0;
    while (!ss.eof()) {
      std::string minor, value;
      ss >> minor;
      if (minor.find('=') != std::string::npos) {
        std::stringstream ss2(minor);
        std::getline(ss2, minor, '=');
        std::getline(ss2, value);
      } else {
        ss >> value;
      }
      //DEBUG("value " << value);
      if (minor.empty()) {
        break;  // skip trailing whitespace
      }
      DEBUG("value " << value);
      if (value[0] == '?') {
        DEBUG("value size " << value.size());
        if (static_cast<int>(value.size()) > 1) {
          value = value.substr(1);
          DEBUG("value " << value);
        }
      }
      ASSERT(!value.empty(), "Error parsing text file on line: \"" << ss.str()
        << "\". Line syntax typically requires an odd number of space-separated "
        << "strings (e.g., Object key0 value0 ... keyN valueN. This error "
        << "typically occurs when one of a key/value pair is missing.");
      DEBUG("major " << major << " minor " << minor << " value " << value);
      if (major == "set_variable") {
        DEBUG("setting variable");
        WARN("Deprecated set_variable->Let");
        ASSERT(variables, "set_variable found when not expected");
        (*variables)[minor] = value;
        *assign_to_list = false;
      } else if (variables && variables->count(value) > 0) {
        DEBUG("using variable");
        args[minor] = (*variables)[value];
      } else {
        DEBUG("no variable: " << value << " sz " << value.size());
        if (value.size() > 7) {
          DEBUG(value.substr(0, 7));
          if (value.substr(0, 7) == "/feasst") {
            DEBUG("replaced: " << value);
            value.replace(0, 7, FEASST_INSTALL_DIR);
          }
        }
        ASSERT(args.find(minor) == args.end(), "All arguments must be unique, but"
          << " the following argument of " << major << " was already used: "
          << minor);
        args[minor] = value;
      }
      ++num_pairs;
      ASSERT(num_pairs < 1e8, "reached maximum number of pairs");
    }
  }
  return std::pair<std::string, argtype>(major, args);
}

argtype line_to_argtype(const std::string line) {
  argtype args;
  if (is_found_in(line, "=")) {
    for (const std::string& pair : split(line, ' ')) {
      const std::vector<std::string> vals = split(pair, '=');
      ASSERT(static_cast<int>(vals.size()) == 2, "TrialGrowFile expects " <<
        "argument pairs to be space-separated, and each pair to be = " <<
        "separated.");
      args[vals[0]] = vals[1];
    }
  } else {
    WARN("Deprecated space separated pairs. Use an equal sign between the " <<
      "name and value of each argument pair.");
    std::stringstream ss(line);
    std::string key, value;
    while (!ss.eof()) {
      ss >> key >> value;
      args[key] = value;
    }
  }
  return args;
}

void replace_value(const std::string search, const std::string replace,
                   arglist * args) {
  for (std::pair<std::string, argtype>& pair1 : *args) {
    for (auto& pair2 : pair1.second) {
      if (pair2.second == search) {
        pair2.second = replace;
      }
    }
  }
}

void replace_in_value(const std::string& from, const std::string& to,
                     arglist * args) {
  for (std::pair<std::string, argtype>& pair1 : *args) {
    for (auto& pair2 : pair1.second) {
      replace(from, to, &pair2.second);
    }
  }
}

std::vector<double> parse_dimensional(const std::string& key, argtype * args,
    const int max) {
  int dim = 0;
  std::stringstream ss;
  ss << key << dim;
  std::vector<double> data;
  while (used(ss.str(), *args)) {
    data.push_back(dble(ss.str(), args));
    ASSERT(dim < max, "dim: " << dim << " is > max: " << max);
    ++dim;
    ss.str("");
    ss << key << dim;
  }
  return data;
}

std::vector<std::pair<std::string, std::vector<std::string> > > parse_for(const argtype& cargs) {
  argtype args = cargs;
  std::vector<std::pair<std::string, std::vector<std::string> > > replace;
  if (args.size() == 0) return replace;
  ASSERT(args.size() == 1, "For takes only one argument at most, but " <<
    args.size() << " arguments were given.");
//  DEBUG(args.begin()->first);
  const std::string vars = args.begin()->first;
  DEBUG("vars " << vars);
  const std::string vals = str(vars, &args);
  DEBUG("vals " << vals);
  ASSERT(!is_found_in(vars, ","), "For loop variable definitions are separated "
    << "by a colon, but a comma was found instead.");
  const std::vector<std::string> nvars = split(vars, ':');
  DEBUG("nvars " << feasst_str(nvars));
  const std::vector<std::string> mvals = split(vals, ',');
  DEBUG("mvals " << feasst_str(mvals));
  replace.resize(nvars.size());
  for (int ivar = 0; ivar < static_cast<int>(nvars.size()); ++ivar) {
    replace[ivar].first = nvars[ivar];
    replace[ivar].second.resize(mvals.size());
    ASSERT(nvars[ivar].front() == '[',
      "For variables should begin with a \"[\" character");
    ASSERT(nvars[ivar].back() == ']',
      "For variables should end with a \"]\" character");
  }
  std::vector<std::vector<std::string> > mnvals;
  for (const std::string& mval : mvals) {
    mnvals.push_back(split(mval, ':'));
    ASSERT(mnvals.back().size() == nvars.size(), "While parsing For arguments, "
      << "the number of given colon-separated variables:" << nvars.size() <<
      " in the line:\"" << vars << "\" did not equal the number "
      << "of given colon-separated values:" << mnvals.back().size() <<
      " in the line:\"" << mval);
  }
  for (int ivar = 0; ivar < static_cast<int>(nvars.size()); ++ivar) {
    for (int ival = 0; ival < static_cast<int>(mvals.size()); ++ival) {
      std::string * mnv = &mnvals[ival][ivar];
      if ((*mnv)[0] == '?' && static_cast<int>(mnv->size()) > 1) {
        *mnv = mnv->substr(1);
      }
      replace[ivar].second[ival] = *mnv;
    }
  }
  for (int ivar = 0; ivar < static_cast<int>(nvars.size()); ++ivar) {
    DEBUG(ivar << " " << replace[ivar].first << " " << feasst_str(replace[ivar].second));
  }
  feasst_check_all_used(args);
  return replace;
}

arglist parse_if(const arglist& list, int * first_end_if, int *last_if, int *last_else) {
  *first_end_if = -1;
  *last_if = -1;
  *last_else = -1;
  const int num = static_cast<int>(list.size());
  bool is_true = false;
  for (int iarg = 0; iarg < num; ++iarg) {
    const std::pair<std::string, argtype>& arg = list[iarg];
    DEBUG("arg.first " << arg.first);
    if (arg.first == "EndIf" || arg.first == "Endif") {
      *first_end_if = iarg;
      break;
    }
    if (arg.first == "Else") {
      ASSERT(*last_else < *last_if, "found an Else without a corresponding If.");
      *last_else = iarg;
    } else if (arg.first == "If") {
      *last_if = iarg;
      DEBUG("true " << str(arg.second));
      const auto pair1 = arg.second.find("defined");
      const auto pair2 = arg.second.find("undefined");
      ASSERT(arg.second.size() == 1 && (pair1 != arg.second.end() || pair2 != arg.second.end()),
       "If syntax is \"If defined=?optional\" or \"If undefined=?opt\". " <<
       "If statement requires only one condition.");
      if (pair1 != arg.second.end()) {
        if (pair1->second == "?") {
          is_true = false;
        } else {
          is_true = true;
        }
      } else if (pair2 != arg.second.end()) {
        if (pair2->second == "?") {
          is_true = true;
        } else {
          is_true = false;
        }
      }
    }
  }
  if (*first_end_if != -1) {
    arglist nwlst;
    ASSERT(*last_if != -1, "Found an EndIf without a corresponding If");
    for (int iarg = 0; iarg < *last_if; ++iarg) {
      nwlst.push_back(list[iarg]);
      DEBUG(iarg << " added: " << nwlst.back().first);
    }
    // process the if substitution
    int start, stop;
    if (is_true) {
      start = *last_if;
      stop = *first_end_if;
      if (*last_else != -1) stop = *last_else;
      DEBUG("start: " << start);
    } else {
      start = *last_else;
      if (*last_else == -1) start = *first_end_if;
      DEBUG("start: " << start);
      stop = *first_end_if;
    }
    for (int iarg = start + 1; iarg < stop; ++iarg) {
      nwlst.push_back(list[iarg]);
      DEBUG(iarg << " added: " << nwlst.back().first);
    }
    for (int iarg = *first_end_if + 1; iarg < num; ++iarg) {
      nwlst.push_back(list[iarg]);
      DEBUG(iarg << " added: " << nwlst.back().first);
    }
    return nwlst;
  } else {
    ASSERT(*last_if == -1, "Found an If without a corresponding EndIf");
    ASSERT(*last_else == -1, "Found an Else without a corresponding EndIf");
    return list;
  }
}

arglist expand_for(const arglist& list, int * first_end_for, int *last_for) {
  std::vector<std::pair<std::string, std::vector<std::string> > > replace;
  *first_end_for = -1;
  *last_for = -1;
  const int num = static_cast<int>(list.size());
  // find an inner most for loop
  for (int iarg = 0; iarg < num; ++iarg) {
    DEBUG(iarg << " " << num);
    const std::pair<std::string, argtype>& arg = list[iarg];
    DEBUG(iarg << " " << arg.first);
    DEBUG(arg.first);
    if (arg.first == "EndFor" || arg.first == "Endfor") {
      *first_end_for = iarg;
      DEBUG("breaking");
      break;
    }
    if (arg.first == "For") {
      *last_for = iarg;
      replace = parse_for(arg.second);
    }
  }
  DEBUG("first " << *first_end_for);
  DEBUG("last " << *last_for);
  if (*first_end_for != -1) {
    // expand the inner most for loop
    arglist nwlst;
    ASSERT(*last_for != -1, "Found an EndFor without a corresponding For");
    // add lines before the last for
    for (int iarg = 0; iarg < *last_for; ++iarg) {
      nwlst.push_back(list[iarg]);
    }
    // process the for loop substitution
    const int copies = static_cast<int>(replace[0].second.size());
    for (int copy = 0; copy < copies; ++copy) {
      for (int iarg = *last_for + 1; iarg < *first_end_for; ++iarg) {
        nwlst.push_back(list[iarg]);
        std::vector<std::string> keys, values;
        for (auto& p : list[iarg].second) {
          keys.push_back(p.first);
          values.push_back(p.second);
        }
        for (const auto& vars : replace) {
          feasst::replace(vars.first, vars.second[copy], &nwlst.back().first);
          for (std::string& str : keys) {
            feasst::replace(vars.first, vars.second[copy], &str);
          }
          for (std::string& str : values) {
            feasst::replace(vars.first, vars.second[copy], &str);
          }
        }
        nwlst.back().second.clear();
        for (int pairs = 0; pairs < static_cast<int>(keys.size()); ++pairs) {
          nwlst.back().second.insert({keys[pairs], values[pairs]});
        }
      }
    }

    // add the rest of the lines
    for (int iarg = *first_end_for + 1; iarg < num; ++iarg) {
      nwlst.push_back(list[iarg]);
    }
    return nwlst;
  } else {
    DEBUG("not expanding For");
    ASSERT(*last_for == -1, "Found a For without a corresponding EndFor");
    return list;
  }
}

std::vector<arglist> parse_mcs(std::istream& is, argtype variables) {
  std::vector<arglist> lists;
  arglist list;
  std::string line;
  while (std::getline(is, line)) {
    if (!line.empty() && line[0] != '#') {
      if (line == "MonteCarlo") {
        lists.push_back(list);
        list = arglist();
      } else {
        bool assign_to_list = true;
        std::pair<std::string, argtype> line_pair = parse_line(line, &variables, &assign_to_list);
        if (!line_pair.first.empty() && assign_to_list) {
          list.push_back(line_pair);
        }
      }
    }
  }
  lists.push_back(list);
  std::vector<arglist> nwlsts;
  for (arglist& list : lists) {
    // parse for
    int first_end_for, last_for;
    arglist nwlst = expand_for(list, &first_end_for, &last_for);
    int attempt = 0;
    while (first_end_for != -1) {
      nwlst = expand_for(nwlst, &first_end_for, &last_for);
      ++attempt;
      ASSERT(attempt < 1e5, "error");
    }
    DEBUG("after parsing for, but before if");
    for (const auto& nw : nwlst) {
      DEBUG(nw.first);
    }

    DEBUG("parsing If");
    int first_end_if, last_if, last_else;
    arglist nwlst2 = parse_if(nwlst, &first_end_if, &last_if, &last_else);
    attempt = 0;
    while (first_end_if != -1) {
      DEBUG("parsing more If");
      nwlst2 = parse_if(nwlst2, &first_end_if, &last_if, &last_else);
      ++attempt;
      ASSERT(attempt < 1e5, "error");
    }
    nwlsts.push_back(nwlst2);
  }
  return nwlsts;
}

}  // namespace feasst
