#!/usr/bin/env python3

from fishing_game_core.game_tree import Node
from fishing_game_core.game_tree import State
from fishing_game_core.player_utils import PlayerController
from fishing_game_core.shared import ACTION_TO_STR

import numpy as np
from math import inf

class PlayerControllerHuman(PlayerController):
    def player_loop(self):
        """
        Function that generates the loop of the game. In each iteration
        the human plays through the keyboard and send
        this to the game through the sender. Then it receives an
        update of the game through receiver, with this it computes the
        next movement.
        :return:
        """

        while True:
            # send message to game that you are ready
            msg = self.receiver()
            if msg["game_over"]:
                return


class PlayerControllerMinimax(PlayerController):

    def __init__(self):
        super(PlayerControllerMinimax, self).__init__()

    def player_loop(self):
        """
        Main loop for the minimax next move search.
        :return:
        """

        # Generate game tree object
        first_msg = self.receiver()
        # Initialize your minimax model
        model = self.initialize_model(initial_data=first_msg)

        while True:
            msg = self.receiver()

            # Create the root node of the game tree
            node = Node(message=msg, player=0)

            # Possible next moves: "stay", "left", "right", "up", "down"
            best_move = self.search_best_next_move(
                model=model, initial_tree_node=node)

            # Execute next action
            self.sender({"action": best_move, "search_time": None})

    def initialize_model(self, initial_data: dict):
        """
        Initialize your minimax model 
        :param initial_data: Game data for initializing minimax model
        :type initial_data: dict
        :return: Minimax model
        :rtype: object
        Sample initial data:
        { 'fish0': {'score': 11, 'type': 3}, 
          'fish1': {'score': 2, 'type': 1}, 
          ...
          'fish5': {'score': -10, 'type': 4},
          'game_over': False }
        Please note that the number of fishes and their types is not fixed between test cases.
        """
        # TODO: EDIT THIS METHOD TO RETURN A MINIMAX MODEL ###

        return None

    def search_best_next_move(self, model: object, initial_tree_node: Node):
        """
        Use your minimax model to find best possible next move for player 0 (green boat)
        :param model: Minimax model
        :type model: object
        :param initial_tree_node: Initial game tree node 
        :type initial_tree_node: game_tree.Node 
            (see the Node class in game_tree.py for more information!)
        :return: either "stay", "left", "right", "up" or "down"
        :rtype: str
        """

        # TODO: EDIT THIS METHOD TO RETURN BEST NEXT POSSIBLE MODE FROM MINIMAX MODEL ###
        
        # NOTE: Don't forget to initialize the children of the current node 
        #       with its compute_and_get_children() method!
        

        #if there is fish on hook, don't search
        if initial_tree_node.state.player_caught[0] != -1:
            return ACTION_TO_STR[1]

        alpha = -inf
        beta = inf
        max_depth = 6
        best_move = None
        best_score = -inf

        '''
        iterativly deepens the search tree
        until max depth is reached
        *Why it's good: 
            in most cases it reduces the number of branches you have to search untill the best move is found
            because the # of braches incraeses exp with the branching factor (aka how many possible S from a given S) 

            the # nodes you revisit (cause look at the same depth several time) is insignifacnt
            especially with alpha beta pruning cause u can prune future branches
        '''
        for i in range(1, max_depth):
            move, score = self.minimax(initial_tree_node, i, 0, alpha, beta)

            alpha = max(alpha, score)

            if score > best_score:
                best_score = score
                best_move = move        

        # 0 is stay, 1 i s up, 2 is down, 
        # 3 is left, 4 is right
        return ACTION_TO_STR[best_move]


    def minimax(self, node: Node, depth: int, player: int, alpha: float, beta: float): # the node we looking at (positions)

        # reaches the max depth or terminal state (no fish left) 
        if depth == 0 or len(node.state.fish_positions) == 0:
            return None, self.heuristic(node.state, player) # retun static evaluation

        #initialize the children of the current node
        children = node.compute_and_get_children()

        # TODO:
        # order the children based on the highest possibility of good score
        # so more pruning can be done
        #ordered = sorted(children, key=lambda n: self.heuristic(n.state, player))
        #ordered = children

        # inspired by https://www.youtube.com/watch?v=l-hh51ncgDI

        if player == 0: # player A - maximizing player

            # bestPossible = −∞
            max_score = -inf 
            best_move = None
            for child in children:

                # v = minimax(child,B)
                # v = max( v, alphabeta(child, depth-1, α, β, B) )
                _, score = self.minimax(child, depth - 1, player^1, alpha, beta) #player XOR 1: switch player each depth
                
                # bestPossible = max(bestPossible, v)      [max_score = max(max_score, score)]
                if score > max_score:
                    max_score = score
                    best_move = child.move

                # α = max(α, v)
                alpha = max(alpha, score)

                #if β ≤ α   break
                if beta <= alpha:
                    break # beta prune

            return best_move, max_score

        else: # player B - minimizing player
                       
            # bestPossible = ∞
            min_score = inf
            best_move = None
            for child in children:

                # v = minimax(child,A)
                # v = min(v , alphabeta(child, depth-1, α, β, A))
                _, score = self.minimax(child, depth - 1, player^1, alpha, beta)

                # bestPossible = min(bestPossible, v)        [min_score = min(min_score, score)]
                if score < min_score:
                    min_score = score
                    best_move = child.move

                # β = min(β, v)
                beta = min(beta, score)

                # i f β ≤ α break 
                if beta <= alpha:
                    break # alpha prune

            return best_move, min_score

    
    ### naive heu function : ν(A,s) = Score(Green-A boat)−Score(Red-B boat) 
    def naive_heuristic(self, state: State) -> float:
        return state.player_scores[0] - state.player_scores[1]


    def heuristic(self, state: State, player: int):
        # initial score is from naive heu
        score = self.naive_heuristic(state)

        good_fish = self.get_good_fish(state.fish_scores, state.fish_positions)

        # If all good fish are taken, only worry about how many points the game will end at (return naive heu)
        # we can't catch any more fish, just return the score from naive heu
        if len(good_fish) == 0:
            return score
        
        bad_fish = self.get_bad_fish(state.fish_scores, state.fish_positions)

        p_hook = state.hook_positions[0] # get the pos of the players hook
        o_hook = state.hook_positions[1]

        # Worse position if close to bad fish
        # if any bad fish is less than 2 unit dist away, decrement/punish the score/player
        if any(self.dist(state.fish_positions[fish], p_hook) < 2 for fish in bad_fish.keys()):
            score -= 1
        
        # Better position of fish hooked
        # if caught a fish, increment/reward score/player
        if state.player_caught[0] != -1:
            score += 1

        # TODO: Base rest of score on distance to good fish. 
        # ### Didnt have time to do 

        return score


    ## HELPER FUNCTIONS ##

    # to get the distance from fish to hook
  
    def dist(self, hook: tuple, fish: tuple) -> float:
        return np.sqrt((hook[0] * fish[0]) + (hook[1] * fish[1]))
   
    '''
    def dist(self, p_hook: tuple, o_hook: tuple, fish: tuple) -> float:
        if p_hook[0] < o_hook[0] and fish[0] > o_hook[0]:
            delta_x = p_hook[0] + 19 - fish[0]
        else:
            delta_x = p_hook[0] - fish[0]
        dist = np.sqrt(np.power(delta_x, 2) + np.power(p_hook[1] - fish[1], 2))
        return max(1, dist)
    '''


    '''
    fishes_scores: map fish key -> fish score
    fishes_pos: map fish key -> fish pos

    for every key k and value v (score) in fishes_score
    if 
        v > 0 (it's not a bad fish)
        and 
        key is in fishes_pos (if fish is on hook, it's in fishes_score but not fishes_pos)
    keep it in dict
    '''
    def get_good_fish(self, fishes_scores: dict, fishes_pos: dict):
        return {k : v for k, v in fishes_scores.items() if v > 0 and k in fishes_pos.keys()}

    #opposite of get_good_fish
    def get_bad_fish(self, fishes_scores: dict, fishes_pos: dict):
        return {k : v for k, v in fishes_scores.items() if v < 0 and k in fishes_pos.keys()}





    '''
    TODO: 
    def get_closest_fish(self, fishes_pos: dict, hook: tuple) -> tuple:
        fish_key = None
        fish_pos = (inf, inf)
        closest_dist = inf
        closest_dist = inf

        for k in fishes_pos.keys():
            dist = self.dist(fishes_pos[k], hook)
            if dist < closest_dist:
                fish_key = k
                fish_pos = fishes_pos[k]
                closest_dist = dist
        return fish_key, fish_pos, closest_dist
    '''