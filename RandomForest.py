# -*- coding: utf-8 -*-
"""
Created on Tue Aug 26 22:12:03 2025

@author: manim
"""
from albero_decisorio import DecisionTree
import numpy as np
from collections import Counter


class RandomForest:
    def __init__(self, n_trees=100, max_depth=5, min_samples_split=2, n_features=None):
        self.n_trees = n_trees
        self.max_depth = max_depth
        self.min_samples_split = min_samples_split
        self.n_features = n_features
        self.trees = []
        
        
    def fit(self, X, y):
        self.trees = []
        for _ in range(self.n_trees):
            tree = DecisionTree(max_depth=self.max_depth, min_samples_split=self.min_samples_split, n_features=self.n_features)
            X_sample, y_sample = self._bootstrap(X, y)
            tree.fit(X_sample, y_sample)
            self.trees.append(tree)
            
    def _bootstrap(self, X, y):
        n_samples = X.shape[0]
        idxs = np.random.choice(n_samples, n_samples, replace=True)
        return X[idxs], y[idxs]
         
    def _most_common_label(self, y):
        counter = Counter(y)
        value = counter.most_common(1)[0][0]
        return value
        
    def predict(self, X):
        predictions = np.array([tree.predict(X) for tree in self.trees])
        tree_predictions = np.swapaxes(predictions, 0, 1)
        # predictions = np.array([self._most_common_label(pred) for pred in tree_predictions])
        predictions = np.array([pred.mean() for pred in tree_predictions])
        stds = np.array([pred.std() for pred in tree_predictions])
        return predictions, stds