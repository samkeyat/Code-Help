# -*- coding: utf-8 -*-
"""
Created on Mon Jul 19 16:03:55 2021

@author: PC
"""

import numpy as np
import pandas as pd
import random

def play_again():
    print("0 : Quit")
    print("1 : Play Game")
    
def play_game():
    random_num = random.randint(0,100)
#    print(random_num)
    guesses = []
    counts = 0

    while True:

        try:
            guess = int(input("Please enter a number from 0 --- 100: "))
            counts +=1
            if guess == random_num:
                print("You Win")
                break

            elif guess < random_num:
                print("Too low")

            else:
                print("Too High")
            guesses.append(int(guess))

        except ValueError:
            print("Please enter a number: ")

    print("Your took", counts, "guesses to win")
    print("your unsuccessful guesses are: ", guesses)
    
    if counts <= 5:
        print("Excellent! You're a genius")
        play_again()

    elif counts <=10:
        print("Good! But you need more work to be a genius")
        play_again()

    else:
        print("You dummy! You're hopeless")
        play_again()

def main_func():

    while True:

        try:
            user_input = int(input())

            if user_input == 0:
                print("Game End")
                break
            elif user_input == 1:
                play_game()
            else:
                print("press 0 or 1")
        except ValueError:
            print("press 0 or 1")
play_again()
main_func()