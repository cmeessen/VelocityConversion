# -*- coding: utf-8 -*-
###############################################################################
#                 Copyright (C) 2017-2022 by Christian Meeßen                 #
#                                                                             #
#                   This file is part of VelocityConversion                   #
#                                                                             #
#  VelocityConversion is free software: you can redistribute it and/or modify #
#    it under the terms of the GNU General Public License as published by     #
#           the Free Software Foundation version 3 of the License.            #
#                                                                             #
#  VelocityConversion is distributed in the hope that it will be useful, but  #
#         WITHOUT ANY WARRANTY; without even the implied warranty of          #
#      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU       #
#                  General Public License for more details.                   #
#                                                                             #
#     You should have received a copy of the GNU General Public License       #
#  along with VelocityConversion. If not, see <http://www.gnu.org/licenses/>. #
###############################################################################
from VelocityConversion import MantleConversion


def main():
    Instance = MantleConversion()
    Instance.ReadArgs()
    Instance.LoadFile()
    Instance.FillTables()
    Instance.CalcPT()
    Instance.SaveFile()


if __name__ == "__main__":
    main()
